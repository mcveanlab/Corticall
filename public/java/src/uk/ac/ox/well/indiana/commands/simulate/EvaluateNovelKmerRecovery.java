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

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 41;

    @Output
    public PrintStream out;

    private class KmerEntry {
        public String variantId;
        public boolean found = false;

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

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                CortexKmer ck = new CortexKmer(seq.substring(i, i + KMER_SIZE));

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
            int start = Integer.valueOf(te.get("start")) + 1 - (KMER_SIZE-1);
            int stop  = Integer.valueOf(te.get("stop")) + 1 + (KMER_SIZE-1);

            if (start < 0) { start = 0; }
            if (stop > REFERENCE.getSequence(chrom).length()) {
                stop = REFERENCE.getSequence(chrom).length();
            }

            String seq = new String(REFERENCE.getSubsequenceAt(chrom, start, stop).getBases());

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                CortexKmer ck = new CortexKmer(seq.substring(i, i + KMER_SIZE));

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

            if (isNovelKmer) {
                novelKmers++;

                //log.info("  novel {}", cr);

                if (variantKmers.containsKey(cr.getCortexKmer())) {
                    variantKmers.get(cr.getCortexKmer()).found = true;
                }
            }
        }
        log.info("  {} novel kmers", novelKmers);

        log.info("Summarizing...");
        Map<String, Integer> expected = new HashMap<String, Integer>();
        Map<String, Integer> actual = new HashMap<String, Integer>();

        for (CortexKmer ck : variantKmers.keySet()) {
            String variantId = variantKmers.get(ck).variantId;
            Map<String, String> te = variants.get(variantId);
            Map<String, String> variantInfo = parseInfoField(te.get("info"));

            String denovo = variantInfo.get("denovo");
            String nahr = variantInfo.get("nahr");

            String vclass = denovo;
            if (denovo.equals("unknown") && !nahr.equals("unknown")) {
                vclass = nahr;
            }

            if (!expected.containsKey(denovo)) {
                expected.put(denovo, 0);
                actual.put(denovo, 0);
            }

            expected.put(denovo, expected.get(denovo) + 1);
            if (variantKmers.get(ck).found) {
                actual.put(denovo, actual.get(denovo) + 1);
            }
        }

        log.info("  exp={}", Joiner.on(", ").withKeyValueSeparator("=").useForNull("null").join(expected));
        log.info("  act={}", Joiner.on(", ").withKeyValueSeparator("=").useForNull("null").join(actual));

        TableWriter tw = new TableWriter(out);

        for (String vclass : expected.keySet()) {
            Map<String, String> te = new LinkedHashMap<String, String>();

            te.put("class", vclass);
            te.put("expected", String.valueOf(expected.get(vclass)));
            te.put("actual", String.valueOf(actual.get(vclass)));

            tw.addEntry(te);
        }
    }
}
