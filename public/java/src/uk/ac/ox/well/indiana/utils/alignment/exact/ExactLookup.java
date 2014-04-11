package uk.ac.ox.well.indiana.utils.alignment.exact;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.util.Interval;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class ExactLookup {
    private Map<String, String> ref = new HashMap<String, String>();
    private Map<Integer, Map<String, Set<Interval>>> lt = new HashMap<Integer, Map<String, Set<Interval>>>();

    private static final int SMALL = (int) Math.pow(2, 8) - 1;
    private static final int MEDIUM = (int) Math.pow(2, 10) - 1;
    private static final int LARGE = (int) Math.pow(2, 12) - 1;

    public ExactLookup(ReferenceSequenceFile sequences) {
        hashSequences(sequences);
    }

    private void hashSequences(ReferenceSequenceFile sequences) {
        lt.put(SMALL, new HashMap<String, Set<Interval>>());
        lt.put(MEDIUM, new HashMap<String, Set<Interval>>());
        lt.put(LARGE, new HashMap<String, Set<Interval>>());

        ReferenceSequence rseq;
        while ((rseq = sequences.nextSequence()) != null) {
            String name = rseq.getName().split("\\s+")[0];
            String seq = new String(rseq.getBases());

            if (!ref.containsKey(name)) {
                ref.put(name, seq);
            }

            System.out.println(name + "...");

            for (int k : new int[]{SMALL, MEDIUM, LARGE}) {
                for (int i = 0; i <= seq.length() - k; i++) {
                    String kmer = seq.substring(i, i + k);

                    Interval interval = new Interval(name, i, i + k);

                    if (!lt.get(k).containsKey(kmer)) {
                        lt.get(k).put(kmer, new HashSet<Interval>());
                    }

                    lt.get(k).get(kmer).add(interval);
                }
            }
        }
    }

    public Interval find(String s) {
        String fw = s.replaceAll("-", "");

        int size;

        if (fw.length() >= SMALL && fw.length() < MEDIUM) {
            size = SMALL;
        } else if (fw.length() >= MEDIUM && fw.length() < LARGE) {
            size = MEDIUM;
        } else {
            size = LARGE;
        }

        if (s.length() >= size) {
            String kmer = s.substring(0, size);

            if (lt.get(size).containsKey(kmer) && lt.get(size).get(kmer).size() == 1) {
                return lt.get(size).get(kmer).iterator().next();
            }
        }

        Set<Interval> homes = new HashSet<Interval>();
        for (String name : ref.keySet()) {
            int index = ref.get(name).indexOf(s);

            if (index >= 0) {
                homes.add(new Interval(name, index, index + s.length()));
            }
        }

        if (homes.size() == 1) {
            return homes.iterator().next();
        }

        return null;
    }
}
