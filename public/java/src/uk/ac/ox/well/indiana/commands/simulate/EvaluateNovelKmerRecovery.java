package uk.ac.ox.well.indiana.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
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

    @Argument(fullName="ref0", shortName="r0", doc="Reference parent 0")
    public FastaSequenceFile REF0;

    @Argument(fullName="ref1", shortName="r1", doc="Reference parent 1")
    public FastaSequenceFile REF1;

    @Argument(fullName="contigMetrics", shortName="m", doc="Contig metrics", required=false)
    public File METRICS;

    @Argument(fullName="bed", shortName="b", doc="Bed file")
    public File BED;

    @Argument(fullName="pedigreeGraph", shortName="g", doc="Pedigree graph")
    public CortexGraph GRAPH;

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

    private Set<CortexKmer> loadNovelKmersFromReferences() {
        Set<CortexKmer> novelKmers = new HashSet<CortexKmer>();

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

        while ((rseq = REF0.nextSequence()) != null) {
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

        while ((rseq = REF1.nextSequence()) != null) {
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

        for (CortexKmer ck : kmerCount.keySet()) {
            if (kmerCount.get(ck) == 1) {
                novelKmers.add(ck);
            }
        }
        return novelKmers;
    }

    private Set<CortexKmer> loadNovelKmersFromContigMetrics() {
        Set<CortexKmer> novelKmers = new HashSet<CortexKmer>();

        TableReader tr = new TableReader(METRICS);

        for (Map<String, String> te : tr) {
            String seq = te.get("seq");
            String kmerOrigin = te.get("kmerOrigin");
            int kmerSize = seq.length() - kmerOrigin.length() + 1;

            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                CortexKmer ck = new CortexKmer(seq.substring(i, i + kmerSize));

                if (kmerOrigin.charAt(i) == '.') {
                    novelKmers.add(ck);
                }
            }
        }

        return novelKmers;
    }

    @Override
    public void execute() {
        log.info("Loading novel kmers...");

        Set<CortexKmer> novelKmersFromReference = loadNovelKmersFromReferences();
        log.info("  reference: {}", novelKmersFromReference.size());

        Set<CortexKmer> novelKmersFromMetrics = loadNovelKmersFromContigMetrics();
        log.info("    metrics: {}", novelKmersFromMetrics.size());

        Map<String, Map<String, String>> variants = new HashMap<String, Map<String, String>>();

        Map<CortexKmer, KmerEntry> kmerRecovery = new HashMap<CortexKmer, KmerEntry>();
        for (CortexKmer ck : novelKmersFromReference) {
            kmerRecovery.put(ck, new KmerEntry("other"));
        }

        log.info("Loading novel kmers from variants...");
        TableReader tr = new TableReader(BED, new String[] {"chrom", "start", "stop", "info"});
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

                if (kmerRecovery.containsKey(ck)) {
                    kmerRecovery.put(ck, new KmerEntry(infoMap.get("id")));
                }
            }
        }

        int novelVariantKmers = 0;
        for (CortexKmer ck : kmerRecovery.keySet()) {
            if (!kmerRecovery.get(ck).variantId.equals("other")) {
                novelVariantKmers++;
            }
        }

        log.info("  {} novel variant kmers", novelVariantKmers);

        log.info("Looking for novel kmers in metrics...");
        int novelKmerCountExpected = 0;
        int novelKmerCountUnexpected = 0;

        for (CortexKmer ck : novelKmersFromMetrics) {
            if (kmerRecovery.containsKey(ck)) {
                kmerRecovery.get(ck).found = true;

                novelKmerCountExpected++;
            } else {
                kmerRecovery.put(ck, new KmerEntry("unexpected"));
                kmerRecovery.get(ck).found = true;

                novelKmerCountUnexpected++;
            }
        }
        log.info("  {} expected novel kmers", novelKmerCountExpected);
        log.info("  {} unexpected novel kmers", novelKmerCountUnexpected);

        log.info("Summarizing...");
        Map<String, Integer> expected = new HashMap<String, Integer>();
        Map<String, Integer> found = new HashMap<String, Integer>();

        int kmersSeen = 0;
        for (CortexKmer ck : kmerRecovery.keySet()) {
            kmersSeen++;

            String variantId = kmerRecovery.get(ck).variantId;
            Map<String, String> te = variantId.equals("other") || variantId.equals("unexpected") ? null : variants.get(variantId);
            Map<String, String> variantInfo = te == null ? null : parseInfoField(te.get("info"));

            String vclass = variantInfo == null ? variantId : variantInfo.get("denovo");

            if (variantInfo != null && !variantInfo.get("nahr").equals("unknown")) {
                vclass = "NAHR";
            }

            if (vclass.equals("unknown")) {
                if (variantInfo != null) {
                    vclass = "inherited_" + variantInfo.get("type");
                }
            }

            if (!expected.containsKey(vclass)) {
                expected.put(vclass, 0);
                found.put(vclass, 0);
            }

            expected.put(vclass, expected.get(vclass) + 1);
            if (kmerRecovery.get(ck).found) {
                found.put(vclass, found.get(vclass) + 1);
            }
        }
        log.info("  {} kmers", kmersSeen);

        TableWriter tw = new TableWriter(out);

        for (String vclass : expected.keySet()) {
            Map<String, String> te = new LinkedHashMap<String, String>();

            String[] pieces = vclass.split("\\.");

            te.put("class", pieces[0]);
            te.put("expected", vclass.equals("unexpected") ? "0" : String.valueOf(expected.get(vclass)));
            te.put("found", String.valueOf(found.get(vclass)));

            tw.addEntry(te);
        }
    }
}
