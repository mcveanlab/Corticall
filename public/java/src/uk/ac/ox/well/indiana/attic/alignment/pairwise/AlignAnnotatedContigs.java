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
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class AlignAnnotatedContigs extends Module {
    @Argument(fullName="metrics", shortName="m", doc="Contig metrics")
    public File METRICS;

    @Argument(fullName="reference", shortName="r", doc="Reference sequence")
    public IndexedFastaSequenceFile REFERENCE;

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

    private void align(String seq, String candidateSeq) {

    }

    @Override
    public void execute() {
        log.info("Loading reference kmer map...");
        Map<String, Set<Interval>> kmerMap = loadReferenceKmers();

        log.info("Computing alignments...");
        TableReader tr = new TableReader(METRICS);
        int numRecords = 0;
        for (Map<String, String> te : tr) {
            if (numRecords % (tr.size() / 10) == 0) {
                log.info("  {}/{} records", numRecords, tr.size());
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

                    align(seq, candidateSeq);

                    log.debug("");
                    log.debug("    loc: {}", candidateInterval);
                    log.debug("    can: {}", candidateSeq);
                    log.debug("    seq: {}", seq);
                    log.debug("");
                }

                log.debug("-----------");
                log.debug("");
            }
        }
        log.info("  {} records", numRecords);
    }
}
