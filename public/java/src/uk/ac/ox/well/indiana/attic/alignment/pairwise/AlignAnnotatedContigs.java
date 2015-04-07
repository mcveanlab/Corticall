package uk.ac.ox.well.indiana.attic.alignment.pairwise;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.StringUtil;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.math.MoreMathUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.indiana.utils.math.MoreMathUtils.max;

public class AlignAnnotatedContigs extends Module {
    @Argument(fullName="metrics", shortName="m", doc="Contig metrics")
    public File METRICS;

    @Argument(fullName="reference", shortName="r", doc="Reference sequence")
    public IndexedFastaSequenceFile REFERENCE;

    @Argument(fullName="delta", shortName="tpd", doc="Transition probability for first gap")
    public Double DELTA = 0.2;

    @Argument(fullName="epsilon", shortName="tpe", doc="Transition probability for remaining in gap")
    public Double EPSILON = 0.2;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Output
    public PrintStream out;

    private Map<String, Set<Interval>> loadReferenceKmers() {
        Map<String, Set<Interval>> kmerMap = new HashMap<String, Set<Interval>>();

        ReferenceSequence rseq;
        while ((rseq = REFERENCE.nextSequence()) != null) {
            String seq = new String(rseq.getBases());
            String seqName = rseq.getName().split("\\s+")[0];

            if (seqName.equals("Pf3D7_07_v3")) {
                log.info("  {}", seqName);

                for (int i = 0; i < seq.length() - KMER_SIZE; i++) {
                    String skFw = seq.substring(i, i + KMER_SIZE);
                    String skRc = SequenceUtils.reverseComplement(skFw);

                    Interval iFw = new Interval(seqName, i, i + KMER_SIZE, false, skFw);
                    Interval iRc = new Interval(seqName, i, i + KMER_SIZE, true, skRc);

                    if (!kmerMap.containsKey(skFw)) { kmerMap.put(skFw, new HashSet<Interval>()); }
                    kmerMap.get(skFw).add(iFw);

                    if (!kmerMap.containsKey(skRc)) { kmerMap.put(skRc, new HashSet<Interval>()); }
                    kmerMap.get(skRc).add(iRc);
                }
            }
        }

        return kmerMap;
    }

    private Set<Interval> getAlignmentRegionCandidates(Map<String, Set<Interval>> kmerMap, String seq) {
        IntervalTreeMap<Interval> candidates = new IntervalTreeMap<Interval>();

        for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
            String sk = seq.substring(i, i + KMER_SIZE);

            if (kmerMap.containsKey(sk)) {
                for (Interval interval : kmerMap.get(sk)) {
                    Interval candidateRegion;
                    if (interval.isPositiveStrand()) {
                        candidateRegion = new Interval(interval.getSequence(),
                                interval.getStart() - i,
                                interval.getStart() - i + seq.length() - 1,
                                false,
                                "candidate");
                    } else {
                        candidateRegion = new Interval(interval.getSequence(),
                                interval.getEnd() + i - seq.length(),
                                interval.getEnd() + i - 1,
                                true,
                                "candidate");
                    }

                    if (candidates.containsOverlapping(candidateRegion)) {
                        for (Interval oi : candidates.getOverlapping(candidateRegion)) {
                            if (oi.isNegativeStrand() == candidateRegion.isNegativeStrand()) {
                                int start = (oi.getStart() < candidateRegion.getStart()) ? oi.getStart() : candidateRegion.getStart();
                                int end = (oi.getEnd() > candidateRegion.getEnd()) ? oi.getEnd() : candidateRegion.getEnd();

                                Interval ni = new Interval(candidateRegion.getSequence(), start, end, candidateRegion.isNegativeStrand(), candidateRegion.getName());

                                candidates.remove(oi);
                                candidates.put(ni, ni);
                            }
                        }
                    } else {
                        candidates.put(candidateRegion, candidateRegion);
                    }
                }
            }
        }

        Set<Interval> candidateIntervals = new TreeSet<Interval>();
        for (Interval interval : candidates.keySet()) {
            Interval candidateInterval = candidates.get(interval);

            Interval c = new Interval(candidateInterval.getSequence(),
                                      candidateInterval.getStart() + 1,
                                      candidateInterval.getEnd() + 1,
                                      candidateInterval.isNegativeStrand(),
                                      candidateInterval.getName());

            candidateIntervals.add(c);
        }

        return candidateIntervals;
    }

    private String getAlignmentRegionCandidateSequence(Interval candidateInterval) {
        int refLength = REFERENCE.getSequence(candidateInterval.getSequence()).length();

        int start = candidateInterval.getStart() < 0 ? 0 : candidateInterval.getStart();
        int end   = candidateInterval.getEnd() >= refLength ? refLength - 1 : candidateInterval.getEnd();

        String candidateSeqFw = new String(REFERENCE.getSubsequenceAt(candidateInterval.getSequence(), start, end).getBases());
        String candidateSeqRc = SequenceUtils.reverseComplement(candidateSeqFw);

        return (candidateInterval.isNegativeStrand()) ? candidateSeqRc : candidateSeqFw;
    }

    private void printMatrix(String query, String target, double[][] m) {
        log.info("   {}", Joiner.on(" ").join(Arrays.copyOfRange(target.split(""), 1, target.length() - 1)));

        for (int i = 0; i < m.length; i++) {
            List<Double> values = new ArrayList<Double>();
            for (int j = 0; j < m[0].length; j++) {
                values.add(m[i][j]);
            }

            log.info("{} {}", query.charAt(i), Joiner.on(" ").join(values));
        }
    }

    private int whichMax(double... values) {
        int maxIndex = -1;
        double maxValue = Double.MIN_VALUE;

        for (int i = 0; i < values.length; i++) {
            if (values[i] > maxValue) {
                maxIndex = i;
                maxValue = values[i];
            }
        }

        return maxIndex;
    }

    private void viterbi(String query, String target, double[][] tm, double[][] em) {
        double[][] vm = new double[query.length()][target.length()];
        double[][] vi = new double[query.length()][target.length()];
        double[][] vd = new double[query.length()][target.length()];

        vm[0][0] = 1;
        for (int i = 0; i < query.length(); i++) {
            for (int j = 0; j < target.length(); j++) {
                if (i != 0 && j != 0) {
                    vm[i][j] = em[0][0] * MoreMathUtils.max(
                                 (1 - 2*DELTA)*vm[i - 1][j - 1],
                                 (1 - EPSILON)*vi[i - 1][j - 1],
                                 (1 - EPSILON)*vd[i - 1][j - 1]
                    );

                    vi[i][j] = em[1][1] * MoreMathUtils.max(
                                 (DELTA)*vm[i - 1][j],
                                 (EPSILON)*vi[i - 1][j]
                    );

                    vd[i][j] = em[2][2] * MoreMathUtils.max(
                                 (DELTA)*vm[i][j - 1],
                                 (EPSILON)*vd[i][j - 1]
                    );
                }
            }
        }

        StringBuilder qa = new StringBuilder();
        StringBuilder ta = new StringBuilder();

        int i = query.length() - 1;
        int j = target.length() - 1;

        while (i >= 0 && j >= 0) {
            switch (whichMax(vm[i][j], vi[i][j], vd[i][j])) {
                case 0:
                    qa.insert(0, query.charAt(i));
                    ta.insert(0, target.charAt(j));
                    i -= 1;
                    j -= 1;
                    break;
                case 1:
                    qa.insert(0, query.charAt(i));
                    ta.insert(0, "-");
                    i -= 1;
                    break;
                case 2:
                    qa.insert(0, "-");
                    ta.insert(0, target.charAt(j));
                    j -= 1;
                    break;
            }
        }

        //printMatrix(query, target, vm);
        //printMatrix(query, target, vi);
        //printMatrix(query, target, vd);

        log.debug("    swt: {}", ta.toString());
        log.debug("    swq: {}", qa.toString());
    }

    private void align(String query, String target) {
        // Define states
        String[] states = new String[] { "M", "I", "D" };

        // Define emissions
        String[] emissions = new String[] { "+1,+1", "+1,+0", "+0,+1" };

        // Define transition matrix
        double[][] tm = new double[3][3];
        tm[0] = new double[] { 1 - 2*DELTA, DELTA,   DELTA   };
        tm[1] = new double[] { 1 - EPSILON, EPSILON, 0       };
        tm[2] = new double[] { 1 - EPSILON, 0,       EPSILON };

        // Define emission matrix
        double[][] em = new double[3][3];
        em[0] = new double[] { 1, 0, 0 };
        em[1] = new double[] { 0, 1, 0 };
        em[2] = new double[] { 0, 0, 1 };

        viterbi(query, target, tm, em);
    }

    @Override
    public void execute() {
        log.info("Loading reference kmer map...");
        Map<String, Set<Interval>> kmerMap = loadReferenceKmers();

        log.info("Computing alignments...");
        TableReader tr = new TableReader(METRICS);
        int numRecords = 0;
        for (Map<String, String> te : tr) {
            if (numRecords % (tr.size() / 5) == 0) {
                log.info("  {}/{} (~{}%) records", numRecords, tr.size(), String.format("%.1f", (100.0f*numRecords) / (float) tr.size()));
            }
            numRecords++;

            String seq = te.get("seq");
            String kmerOrigin = te.get("kmerOrigin");

            if (te.get("canonicalLocus").contains("Pf3D7_07_v3")) {
                log.debug("  {}: {} {} {}", te.get("contigName"), te.get("canonicalLocus"), te.get("cigarCanonical"), seq.length());
                log.debug("    seq: {}", seq);
                log.debug("     ko: {}", kmerOrigin);

                Set<Interval> candidates = getAlignmentRegionCandidates(kmerMap, seq);

                for (Interval candidateInterval : candidates) {
                    String candidateSeq = getAlignmentRegionCandidateSequence(candidateInterval);

                    log.debug("");
                    log.debug("    loc: {}", candidateInterval);
                    log.debug("    can: {}", candidateSeq);
                    log.debug("    seq: {}", seq);

                    if (seq.length() < 200) {
                        align(seq, candidateSeq);
                    }

                    log.debug("");
                }

                log.debug("-----------");
                log.debug("");
            }
        }
    }
}
