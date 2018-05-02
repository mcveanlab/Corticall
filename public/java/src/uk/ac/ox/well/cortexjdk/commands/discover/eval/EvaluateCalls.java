package uk.ac.ox.well.cortexjdk.commands.discover.eval;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.alignment.sw.SmithWaterman;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.containers.ContainerUtils;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableWriter;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;

import java.io.File;
import java.util.*;

public class EvaluateCalls extends Module {
    //@Argument(fullName="graph", shortName="g", doc="Graph")
    //public CortexGraph GRAPH;

    @Argument(fullName="known", shortName="k", doc="VCF of known variants")
    public ArrayList<VCFFileReader> VCF_KNOWNS;

    @Argument(fullName="novel", shortName="n", doc="VCF of novel variants")
    public ArrayList<VCFFileReader> VCF_NOVELS;

    @Argument(fullName="acct", shortName="a", doc="Novel variants accounting")
    public ArrayList<File> ACCT;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Output
    public File out;

    private class StatsEntry {
        private int numEvents = 0;
        private int someEvidence = 0;
        private int numCorrectBackground = 0;
        private int numCorrectChromosome = 0;
        private int numCorrectStart = 0;
        private int numCorrectStop = 0;
        private int numCorrectType = 0;
        private int numComplete = 0;
        private int numNovelsExpected = 0;
        private int numNovelsFound = 0;
        private int numFalsePositives = 0;

        @Override
        public String toString() {
            return "StatsEntry{" +
                    "numEvents=" + numEvents +
                    ", someEvidence=" + someEvidence +
                    ", numCorrectBackground=" + numCorrectBackground +
                    ", numCorrectChromosome=" + numCorrectChromosome +
                    ", numCorrectStart=" + numCorrectStart +
                    ", numCorrectStop=" + numCorrectStop +
                    ", numCorrectType=" + numCorrectType +
                    ", numComplete=" + numComplete +
                    ", numNovelsExpected=" + numNovelsExpected +
                    ", numNovelsFound=" + numNovelsFound +
                    ", numFalsePositives=" + numFalsePositives +
                    '}';
        }
    }

    private void addEvent(Map<String, StatsEntry> stats, String type) {
        StatsEntry se = stats.containsKey(type) ? stats.get(type) : new StatsEntry();

        se.numEvents++;

        stats.put(type, se);
    }

    private void addEntry(Map<String, StatsEntry> stats, String type, boolean correctBackground, boolean correctChromosome, boolean correctStart, boolean correctStop, boolean correctType, boolean complete, int numNovelsExpected, int numNovelsFound) {
        StatsEntry se = stats.containsKey(type) ? stats.get(type) : new StatsEntry();

        se.someEvidence++;
        se.numCorrectBackground += correctBackground ? 1 : 0;
        se.numCorrectChromosome += correctChromosome ? 1 : 0;
        se.numCorrectStart += correctStart ? 1 : 0;
        se.numCorrectStop += correctStop ? 1 : 0;
        se.numCorrectType += correctType ? 1 : 0;
        se.numComplete += complete ? 1 : 0;
        se.numNovelsExpected += numNovelsExpected;
        se.numNovelsFound += numNovelsFound;

        stats.put(type, se);
    }

    private void addFp(Map<String, StatsEntry> stats, String type) {
        StatsEntry se = stats.containsKey(type) ? stats.get(type) : new StatsEntry();

        se.numFalsePositives++;

        stats.put(type, se);
    }

    @Override
    public void execute() {
        Map<String, StatsEntry> stats = new HashMap<>();

        for (int q = 0; q < VCF_KNOWNS.size(); q++) {
            VCFFileReader vcfKnown = VCF_KNOWNS.get(q);
            VCFFileReader vcfNovel = VCF_NOVELS.get(q);

            TableReader tr = new TableReader(ACCT.get(q), "novel", "ccid");
            Map<CanonicalKmer, String> ids = new HashMap<>();
            for (Map<String, String> te : tr) {
                CanonicalKmer ck = new CanonicalKmer(te.get("novel"));
                String id = te.get("ccid");

                ids.put(ck, id);
            }

            Map<CanonicalKmer, Set<VariantContext>> kmerMap = new HashMap<>();
            Map<VariantContext, Set<VariantContext>> recovered = new LinkedHashMap<>();
            Set<VariantContext> unused = new HashSet<>();

            int rec = 0, unu = 0, nov = 0;
            for (VariantContext kvc : vcfKnown) {
                addEvent(stats, kvc.getAttributeAsString("SIM_TYPE", "unknown"));

                recovered.put(kvc, new HashSet<>());

                String newHap = kvc.getAttributeAsString("NEW_HAP", "");
                for (int i = 0; i <= newHap.length() - KMER_SIZE; i++) {
                    String sk = newHap.substring(i, i + KMER_SIZE);
                    CanonicalKmer ck = new CanonicalKmer(sk);

                    if (!kmerMap.containsKey(ck)) {
                        kmerMap.put(ck, new HashSet<>());
                    }

                    kmerMap.get(ck).add(kvc);
                }
            }

            for (VariantContext nvc : vcfNovel) {
                if (nvc.getAttributeAsInt("flankMappingQuality", 0) == 60 && nvc.getAttributeAsList("alt_loci").size() <= 2) {
                    nov++;

                    Map<VariantContext, Integer> knownVcs = new HashMap<>();

                    String childHap = nvc.getAttributeAsString("CHILD_HAP", "");
                    for (int i = 0; i <= childHap.length() - KMER_SIZE; i++) {
                        String sk = childHap.substring(i, i + KMER_SIZE);
                        CanonicalKmer ck = new CanonicalKmer(sk);

                        if (kmerMap.containsKey(ck)) {
                            for (VariantContext vc : kmerMap.get(ck)) {
                                ContainerUtils.increment(knownVcs, vc);
                            }
                        }
                    }

                    VariantContext bestVc = null;
                    int bestCount = 0;
                    for (VariantContext knownVc : knownVcs.keySet()) {
                        if (knownVcs.get(knownVc) > bestCount) {
                            bestVc = knownVc;
                            bestCount = knownVcs.get(knownVc);
                        }
                    }

                    if (bestVc != null) {
                        recovered.get(bestVc).add(nvc);
                    } else {
                        String ccid = "CC" + nvc.getAttributeAsString("CALL_ID", "");
                        for (CanonicalKmer ck : ids.keySet()) {
                            if (ids.get(ck).equals(ccid)) {
                                log.info("{} {} {}", ccid, ck, nvc);
                            }
                        }

                        unused.add(nvc);
                    }
                }
            }

            for (VariantContext bestVc : recovered.keySet()) {
                log.info("known: {} {}", bestVc.getAttributeAsString("SIM_TYPE", "unknown"), bestVc);

                // background, start_correct, stop_correct, alleles
                //boolean someEvidence = true;
                boolean correctBackground = false;
                boolean correctChromosome = false;
                boolean correctStart = false;
                boolean correctStop = false;
                boolean correctType = false;
                boolean complete = recovered.get(bestVc).size() == 1;

                List<VariantContext> vs = new ArrayList<>(recovered.get(bestVc));
                vs.sort(Comparator.comparingInt(VariantContext::getStart));

                int numNovelsExpected = 0;
                int numNovelsFound = 0;
                for (VariantContext nvc : vs) {
                    log.info("call : {} {}", nvc.getAttributeAsString("SVTYPE", "unknown"), nvc);

                    if (bestVc.getContig().contains(nvc.getAttributeAsString("BACKGROUND", "unknown"))) {
                        correctBackground = true;
                    }

                    if (bestVc.getContig().equals(nvc.getContig())) {
                        correctChromosome = true;
                    }

                    if (Math.abs(bestVc.getStart() - nvc.getStart()) <= 15) {
                        correctStart = true;
                    }

                    if (Math.abs(bestVc.getEnd() - nvc.getEnd()) <= 15) {
                        correctStop = true;
                    }

                    if (bestVc.getType().equals(nvc.getType())) {
                        correctType = true;
                    }

                    String newHap = bestVc.getAttributeAsString("NEW_HAP", "");
                    String childHap = nvc.getAttributeAsString("CHILD_HAP", "");

                    for (int i = 0; i <= newHap.length() - KMER_SIZE; i++) {
                        CanonicalKmer ck = new CanonicalKmer(newHap.substring(i, i + KMER_SIZE));
                        if (ids.containsKey(ck)) {
                            numNovelsExpected++;
                        }
                    }

                    for (int i = 0; i <= childHap.length() - KMER_SIZE; i++) {
                        CanonicalKmer ck = new CanonicalKmer(childHap.substring(i, i + KMER_SIZE));
                        if (ids.containsKey(ck) && !ids.get(ck).equals("absent")) {
                            numNovelsFound++;
                        }
                    }

                    //SmithWaterman sw = new SmithWaterman();
                    //String[] a = sw.getAlignment(childHap, newHap);

                    //log.info("  {}", a[0]);
                    //log.info("  {}", a[1]);
                }

                addEntry(stats, bestVc.getAttributeAsString("SIM_TYPE", "unknown"), correctBackground, correctChromosome, correctStart, correctStop, correctType, complete, numNovelsExpected, numNovelsFound);

                log.info("       back={} chr={} start={} stop={} type={} complete={} numNovelsExpected={} numNovelsFound={}", correctBackground, correctChromosome, correctStart, correctStop, correctType, complete, numNovelsExpected, numNovelsFound);
                log.info("");
            }

            for (VariantContext vcu : unused) {
                //log.info("unused: {}", unused.size());

                String type = "UNKNOWN";
                if (vcu.isSNP()) {
                    type = "SNV";
                } else if (vcu.isSimpleInsertion()) {
                    if (vcu.getAltAlleleWithHighestAlleleCount().length() >= 500) {
                        type = "LARGE_INS";
                    } else {
                        type = "SMALL_INS";
                    }
                } else if (vcu.isSimpleDeletion()) {
                    if (vcu.getReference().length() >= 500) {
                        type = "LARGE_DEL";
                    } else {
                        type = "SMALL_DEL";
                    }
                } else if (vcu.isMNP()) {
                    if (vcu.getReference().length() >= 500) {
                        type = "LARGE_MNP";
                    } else {
                        type = "SMALL_MNP";
                    }
                } else if (vcu.isSymbolicOrSV()) {
                    type = "BREAKPOINT";
                }

                addFp(stats, type);
            }
        }

        for (String type : stats.keySet()) {
            float tp0 = stats.get(type).someEvidence; // tp
            float fn0 = stats.get(type).numEvents - stats.get(type).someEvidence; // fn
            float fp0 = stats.get(type).numFalsePositives; // fp
            float prec0 = tp0 / (tp0 + fp0);
            float rec0 = tp0 / (tp0 + fn0);
            float f10 = 2*(rec0*prec0) / (rec0 + prec0);

            float tp1 = stats.get(type).numComplete; // tp
            float fn1 = stats.get(type).numEvents - stats.get(type).numComplete; // fn
            float fp1 = stats.get(type).numFalsePositives; // fp
            float prec1 = tp1 / (tp1 + fp1);
            float rec1 = tp1 / (tp1 + fn1);
            float f11 = 2*(rec1*prec1) / (rec1 + prec1);

            if (!type.equals("NAHR-DEL")) {
                log.info("{} {} {} {} [tp0={} fn0={} fp0={}] [tp1={} fn1={} fp1={}], f1={} ({})",
                        type,
                        stats.get(type).numNovelsExpected,
                        stats.get(type).numNovelsFound,
                        stats.get(type).numEvents,

                        tp0,
                        fn0,
                        fp0,

                        tp1,
                        fn1,
                        fp1,

                        String.format("%.2f%%", 100.0f*f10),
                        String.format("%.2f%%", 100.0f*f11)
                );

                /*
                log.info("{} {}", type, stats.get(type));
                log.info("{} {} {} {} complete[tp={} fn={} fp={} prec={} rec={} f1={}]",
                        type,
                        stats.get(type).numNovelsExpected,
                        stats.get(type).numNovelsFound,
                        stats.get(type).numEvents,

                        tp0,
                        fn0,
                        fp0,

                        prec0,
                        rec0,
                        f10
                );
                */
            }
        }

    }
}
