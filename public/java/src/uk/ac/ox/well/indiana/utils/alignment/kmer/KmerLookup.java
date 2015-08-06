package uk.ac.ox.well.indiana.utils.alignment.kmer;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CigarUtil;
import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.util.*;

public class KmerLookup {
    private int kmerSize;
    private KmerIndex ki;

    public KmerLookup(File ref) {
        initialize(ref, 47);
    }

    public KmerLookup(File ref, int kmerSize) {
        initialize(ref, kmerSize);
    }

    private void initialize(File ref, int kmerSize) {
        this.kmerSize = kmerSize;

        this.ki = loadIndex(ref, kmerSize);
    }

    private KmerIndex loadIndex(File ref, int kmerSize) {
        return new KmerIndex(ref, kmerSize);
    }

    public KmerIndex.KmerIndexRecord findKmer(String sk) {
        return ki.lookup(sk);
    }

    public List<Set<Interval>> findKmers(String s) {
        List<Set<Interval>> allIntervals = new ArrayList<Set<Interval>>();

        for (int i = 0; i <= s.length() - kmerSize; i++) {
            String sk = s.substring(i, i + kmerSize);

            KmerIndex.KmerIndexRecord kir = findKmer(sk);

            if (kir == null) {
                allIntervals.add(new HashSet<Interval>());
            } else {
                allIntervals.add(kir.getLocations());
            }
        }

        return allIntervals;
    }

    private List<Set<Interval>> combineIntervals(List<Set<Interval>> allIntervals, boolean isNegative) {
        List<Set<Interval>> combinedIntervals = new ArrayList<Set<Interval>>();

        Interval currentInterval = null;
        for (Set<Interval> intervals : allIntervals) {
            if (intervals.size() == 1) {
                Interval interval = intervals.iterator().next();

                if (currentInterval == null) {
                    currentInterval = new Interval(interval.getSequence(), interval.getStart(), interval.getEnd(), isNegative, "none");
                } else {
                    if (currentInterval.getSequence().equals(interval.getSequence()) && interval.abuts(currentInterval)) {
                        currentInterval = new Interval(
                                currentInterval.getSequence(),
                                currentInterval.getStart() < interval.getStart() ? currentInterval.getStart() : interval.getStart(),
                                currentInterval.getEnd()   > interval.getEnd()   ? currentInterval.getEnd()   : interval.getEnd(),
                                isNegative,
                                "."
                        );
                    } else {
                        Set<Interval> combined = new HashSet<Interval>();
                        combined.add(currentInterval);

                        combinedIntervals.add(combined);

                        currentInterval = interval;
                    }
                }
            } else {
                if (currentInterval != null) {
                    Set<Interval> combined = new HashSet<Interval>();
                    combined.add(currentInterval);

                    combinedIntervals.add(combined);

                    currentInterval = null;
                }
            }
        }

        if (currentInterval != null) {
            Set<Interval> combined = new HashSet<Interval>();
            combined.add(currentInterval);

            combinedIntervals.add(combined);
        }

        return combinedIntervals;
    }

    private int computeNumKmersAlignedUniquely(List<Set<Interval>> allIntervals) {
        int numKmersAlignedUniquely = 0;

        for (Set<Interval> intervals : allIntervals) {
            if (intervals.size() == 1) {
                numKmersAlignedUniquely++;
            }
        }

        return numKmersAlignedUniquely;
    }

    public List<Set<Interval>> align(String sFw) {
        String sRc = SequenceUtils.reverseComplement(sFw);

        List<Set<Interval>> allIntervalsFw = findKmers(sFw);
        List<Set<Interval>> allIntervalsRc = findKmers(sRc);
        int scoreFw = computeNumKmersAlignedUniquely(allIntervalsFw);
        int scoreRc = computeNumKmersAlignedUniquely(allIntervalsRc);

        List<Set<Interval>> combinedIntervals = (scoreFw > scoreRc) ? combineIntervals(allIntervalsFw, false) : combineIntervals(allIntervalsRc, true);

        return combinedIntervals;
    }
}
