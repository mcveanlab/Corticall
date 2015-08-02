package uk.ac.ox.well.indiana.utils.alignment.kmer;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;

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

    public String findSimplified(String s) {
        List<Set<Interval>> allIntervals = findKmers(s);
        StringBuilder simplified = new StringBuilder();

        /*
        Map<String, Integer> chrCounts = new HashMap<String, Integer>();
        for (Set<Interval> intervals : allIntervals) {
            Map<String, Integer> cc = new HashMap<String, Integer>();

            for (Interval interval : intervals) {
                ContainerUtils.increment(cc, interval.getSequence());
            }

            ContainerUtils.increment(chrCounts, ContainerUtils.mostCommonKey(cc));
        }
        */

        for (Set<Interval> intervals : allIntervals) {
            simplified.append(intervals.size() == 0 ? "0" : "1");
        }

        return simplified.toString();
    }
}
