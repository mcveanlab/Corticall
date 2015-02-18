package uk.ac.ox.well.indiana.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class EvaluateNovelKmerRecovery extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference")
    public IndexedFastaSequenceFile REFERENCE;

    @Argument(fullName="bed", shortName="b", doc="Bed file")
    public File BED;

    @Argument(fullName="pedigreeGraph", shortName="g", doc="Pedigree graph")
    public CortexGraph GRAPH;

    @Output
    public PrintStream out;

    private class KmerEntry {
        public String variantId;
        public boolean found = false;
        public boolean novel = false;
        public int coverage = 0;

        public KmerEntry(String variantId) {
            this.variantId = variantId;
        }
    }

    private Map<String, String> parseInfoField(String info) {
        Map<String, String> infoMap = new HashMap<String, String>();
        String[] entries = info.split(";");

        for (String entry : entries) {
            String[] fields = entry.split("=");

            infoMap.put(fields[0], fields[1]);
        }

        return infoMap;
    }

    @Override
    public void execute() {
        log.info("Loading reference kmer counts...");
        Map<CortexKmer, Integer> kmerCount = new HashMap<CortexKmer, Integer>();
        ReferenceSequence rseq;
        while ((rseq = REFERENCE.nextSequence()) != null) {
            String seq = new String(rseq.getBases());

            for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                CortexKmer ck = new CortexKmer(seq.substring(i, i + GRAPH.getKmerSize()));

                if (!kmerCount.containsKey(ck)) {
                    kmerCount.put(ck, 1);
                } else {
                    kmerCount.put(ck, kmerCount.get(ck) + 1);
                }
            }
        }
        log.info("  {} kmers", kmerCount.size());

        TableReader tr = new TableReader(BED, new String[] {"chrom", "start", "stop", "info"});

        Map<String, Map<String, String>> variants = new HashMap<String, Map<String, String>>();
        Map<CortexKmer, KmerEntry> variantKmers = new HashMap<CortexKmer, KmerEntry>();

        log.info("Loading kmers from variants...");
        int kmers = 0;
        for (Map<String, String> te : tr) {
            Map<String, String> infoMap = parseInfoField(te.get("info"));

            variants.put(infoMap.get("id"), te);

            String chrom = te.get("chrom");
            int start = Integer.valueOf(te.get("start")) + 1 - (GRAPH.getKmerSize()-1);
            int stop  = Integer.valueOf(te.get("stop")) + 1 + (GRAPH.getKmerSize()-1);

            if (start < 0) { start = 0; }
            if (stop > REFERENCE.getSequence(chrom).length()) {
                stop = REFERENCE.getSequence(chrom).length();
            }

            String seq = new String(REFERENCE.getSubsequenceAt(chrom, start, stop).getBases());

            for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                CortexKmer ck = new CortexKmer(seq.substring(i, i + GRAPH.getKmerSize()));

                if (kmerCount.containsKey(ck) && kmerCount.get(ck) == 1) {
                    variantKmers.put(ck, new KmerEntry(infoMap.get("id")));
                    kmers++;
                }
            }
        }
        log.info("  {} kmers", kmers);

        log.info("Looking for novel kmers in graph...");
        int novelKmers = 0;

        int childColor = 0;
        for (int c = 0; c < GRAPH.getNumColors(); c++) {
            if (GRAPH.getColor(c).getSampleName().contains("child")) {
                childColor = c;
            }
        }

        for (CortexRecord cr : GRAPH) {
            boolean isNovelKmer = cr.getCoverage(childColor) > 0;

            for (int c = 0; c < GRAPH.getNumColors(); c++) {
                if (c != childColor) {
                    isNovelKmer &= cr.getCoverage(c) == 0;
                }
            }

            if (variantKmers.containsKey(cr.getCortexKmer())) {
                variantKmers.get(cr.getCortexKmer()).found = true;
                variantKmers.get(cr.getCortexKmer()).coverage = cr.getCoverage(childColor);
            } else {
                variantKmers.put(cr.getCortexKmer(), new KmerEntry("other"));
                variantKmers.get(cr.getCortexKmer()).found = true;
                variantKmers.get(cr.getCortexKmer()).coverage = cr.getCoverage(childColor);
            }

            if (isNovelKmer) {
                novelKmers++;

                if (variantKmers.containsKey(cr.getCortexKmer())) {
                    variantKmers.get(cr.getCortexKmer()).novel = true;
                }
            }
        }
        log.info("  {} novel kmers", novelKmers);

        log.info("Summarizing...");
        Map<String, Integer> expected = new HashMap<String, Integer>();
        Map<String, Integer> found = new HashMap<String, Integer>();
        Map<String, Integer> novel = new HashMap<String, Integer>();
        //Map<String, Integer> cov5 = new HashMap<String, Integer>();
        //Map<String, Integer> cov10 = new HashMap<String, Integer>();
        //Map<String, Integer> cov15 = new HashMap<String, Integer>();
        //Map<String, Integer> cov20 = new HashMap<String, Integer>();
        //Map<String, Integer> cov25 = new HashMap<String, Integer>();
        //Map<String, Integer> cov30 = new HashMap<String, Integer>();

        for (CortexKmer ck : variantKmers.keySet()) {
            String variantId = variantKmers.get(ck).variantId;
            Map<String, String> te = variantId.equals("other") ? null : variants.get(variantId);
            Map<String, String> variantInfo = variantId.equals("other") ? null : parseInfoField(te.get("info"));

            //String denovo = variantInfo.get("denovo");
            //int length = 0;
            //if (denovo.equals("INS") || denovo.equals("TD") || denovo.equals("STR_EXP")) { length = variantInfo.get("alt").length() - 1; }
            //if (denovo.equals("DEL") || denovo.equals("STR_CON")) { length = variantInfo.get("ref").length() - 1; }
            //if (denovo.equals("INV")) { length = variantInfo.get("ref").length(); }
            //String vclass = denovo + "." + length;

            String vclass = variantId.equals("other") ? "other" : variantInfo.get("denovo");

            if (!expected.containsKey(vclass)) {
                expected.put(vclass, 0);
                found.put(vclass, 0);
                novel.put(vclass, 0);
                //cov5.put(vclass, 0);
                //cov10.put(vclass, 0);
                //cov15.put(vclass, 0);
                //cov20.put(vclass, 0);
                //cov25.put(vclass, 0);
                //cov30.put(vclass, 0);
            }

            expected.put(vclass, expected.get(vclass) + 1);
            if (variantKmers.get(ck).found) {
                found.put(vclass, found.get(vclass) + 1);

                if (variantKmers.get(ck).novel) {
                    novel.put(vclass, novel.get(vclass) + 1);

                    //if (variantKmers.get(ck).coverage >= 5)  { cov5.put(vclass, cov5.get(vclass) + 1); }
                    //if (variantKmers.get(ck).coverage >= 10) { cov10.put(vclass, cov10.get(vclass) + 1); }
                    //if (variantKmers.get(ck).coverage >= 15) { cov15.put(vclass, cov15.get(vclass) + 1); }
                    //if (variantKmers.get(ck).coverage >= 20) { cov20.put(vclass, cov20.get(vclass) + 1); }
                    //if (variantKmers.get(ck).coverage >= 25) { cov25.put(vclass, cov25.get(vclass) + 1); }
                    //if (variantKmers.get(ck).coverage >= 30) { cov30.put(vclass, cov30.get(vclass) + 1); }
                }
            }
        }

        //log.info("    exp={}", Joiner.on(", ").withKeyValueSeparator("=").useForNull("null").join(expected));
        //log.info("  found={}", Joiner.on(", ").withKeyValueSeparator("=").useForNull("null").join(found));
        //log.info("  novel={}", Joiner.on(", ").withKeyValueSeparator("=").useForNull("null").join(novel));

        TableWriter tw = new TableWriter(out);

        for (String vclass : expected.keySet()) {
            Map<String, String> te = new LinkedHashMap<String, String>();

            String[] pieces = vclass.split("\\.");

            te.put("class", pieces[0]);
            te.put("length", (pieces.length == 2) ? pieces[1] : "0");
            te.put("expected", String.valueOf(expected.get(vclass)));
            te.put("found", String.valueOf(found.get(vclass)));
            te.put("novel", String.valueOf(novel.get(vclass)));
            //te.put("cov5", String.valueOf(cov5.get(vclass)));
            //te.put("cov10", String.valueOf(cov10.get(vclass)));
            //te.put("cov15", String.valueOf(cov15.get(vclass)));
            //te.put("cov20", String.valueOf(cov20.get(vclass)));
            //te.put("cov25", String.valueOf(cov25.get(vclass)));
            //te.put("cov30", String.valueOf(cov30.get(vclass)));

            tw.addEntry(te);
        }
    }
}
