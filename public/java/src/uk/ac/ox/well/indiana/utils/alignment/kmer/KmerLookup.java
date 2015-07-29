package uk.ac.ox.well.indiana.utils.alignment.kmer;

import htsjdk.samtools.util.Interval;

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

    public List<Set<Interval>> find(String s) {
        List<Set<Interval>> allIntervals = new ArrayList<Set<Interval>>();

        for (int i = 0; i <= s.length() - kmerSize; i++) {
            String sk = s.substring(i, i + kmerSize);

            KmerIndex.KmerIndexRecord kir = ki.lookup(sk);

            if (kir == null) {
                allIntervals.add(new HashSet<Interval>());
            } else {
                allIntervals.add(kir.getLocations());
            }
        }

        return allIntervals;
    }
}
