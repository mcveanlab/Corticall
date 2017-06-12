package uk.ac.ox.well.indiana.utils.alignment.kmer;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import org.mapdb.BTreeMap;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.Serializer;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexBinaryKmer;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

//import org.mapdb.Fun;

public class KmerLookup {
    private int kmerSize;
    private BTreeMap<long[], int[]> kmerIndex;
    private IndexedFastaSequenceFile ref;

    public KmerLookup(File refFile) {
        initialize(refFile, 47);
    }

    public KmerLookup(File refFile, int kmerSize) {
        initialize(refFile, kmerSize);
    }

    private void initialize(File refFile, int kmerSize) {
        this.kmerSize = kmerSize;
        try {
            this.ref = new IndexedFastaSequenceFile(refFile);

            if (this.ref.getSequenceDictionary() == null) {
                throw new IndianaException("Reference must have a sequence dictionary");
            }
        } catch (FileNotFoundException e) {
            throw new IndianaException("File not found", e);
        }

        loadIndex(refFile);
    }

    private void loadIndex(File ref) {
        File dbFile = new File(ref.getAbsoluteFile() + ".kmerdb");

        if (dbFile.exists()) {
            DB db = DBMaker.fileDB(dbFile)
                    .fileMmapEnable()
                    .closeOnJvmShutdown()
                    .readOnly()
                    .make();

            kmerIndex = db.treeMap("index")
                    .keySerializer(Serializer.LONG_ARRAY)
                    .valueSerializer(Serializer.INT_ARRAY)
                    .counterEnable()
                    .open();
        } else {
            System.err.println("No index for '" + dbFile + "'");
        }
    }

    public String findKmer(Interval interval) {
        ReferenceSequence rseq = ref.getSubsequenceAt(interval.getContig(), interval.getStart(), interval.getEnd());

        if (rseq != null) {
            if (interval.isPositiveStrand()) {
                return rseq.getBaseString();
            } else {
                return SequenceUtils.reverseComplement(rseq.getBaseString());
            }
        }

        return null;
    }

    public Set<Interval> findKmer(String sk) {
        Set<Interval> intervals = new HashSet<>();

        int[] l = kmerIndex.get(new CortexBinaryKmer(sk.getBytes()).getBinaryKmer());

        if (l != null) {
            for (int i = 0; i < l.length - 1; i += 2) {
                String chr = ref.getSequenceDictionary().getSequence(l[i]).getSequenceName();
                int pos = l[i + 1];

                String fw = ref.getSubsequenceAt(chr, pos + 1, pos + sk.length()).getBaseString();
                String rc = SequenceUtils.reverseComplement(fw);

                if (sk.equals(fw)) {
                    intervals.add(new Interval(chr, pos + 1, pos + sk.length()));
                } else if (sk.equals(rc)) {
                    intervals.add(new Interval(chr, pos + 1, pos + sk.length(), true, null));
                }
            }
        }

        return intervals;
    }

    public List<Set<Interval>> findKmers(String s) {
        List<Set<Interval>> allIntervals = new ArrayList<>();

        for (int i = 0; i <= s.length() - kmerSize; i++) {
            String sk = s.substring(i, i + kmerSize);

            allIntervals.add(findKmer(sk));
        }

        return allIntervals;
    }

    private List<Set<Interval>> combineIntervals(List<Set<Interval>> allIntervals, boolean isNegative) {
        List<Set<Interval>> combinedIntervals = new ArrayList<>();

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
                        Set<Interval> combined = new HashSet<>();
                        combined.add(currentInterval);

                        combinedIntervals.add(combined);

                        currentInterval = interval;
                    }
                }
            } else {
                if (currentInterval != null) {
                    Set<Interval> combined = new HashSet<>();
                    combined.add(currentInterval);

                    combinedIntervals.add(combined);

                    currentInterval = null;
                }
            }
        }

        if (currentInterval != null) {
            Set<Interval> combined = new HashSet<>();
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

    private List<Interval> flipStrands(List<Interval> lis) {
        List<Interval> newLis = new ArrayList<>();

        for (Interval li : lis) {
            newLis.add(new Interval(li.getSequence(), li.getStart(), li.getEnd(), true, li.getName()));
        }

        return newLis;
    }

    private Interval closestUniqueAlignment(List<List<Interval>> ka, int index) {
        Interval k1 = null, k2 = null;
        int i1 = Integer.MAX_VALUE, i2 = Integer.MAX_VALUE;

        for (int i = index - 1; i >= 0; i--) {
            if (ka.get(i).size() == 1) {
                k1 = ka.get(i).get(0);
                i1 = i;
            }
        }

        for (int i = index + 1; i < ka.size(); i++) {
            if (ka.get(i).size() == 1) {
                k2 = ka.get(i).get(0);
                i2 = i;
            }
        }

        if      (i1 < i2) { return k1; }
        else if (i2 < i1) { return k2; }

        return null;
    }

    private Interval closestMatchingAlignment(Interval it, List<Interval> ka) {
        if (it != null) {
            Map<Integer, Interval> cands = new TreeMap<>();

            for (Interval ita : ka) {
                if (it.getSequence().equals(ita.getSequence())) {
                    cands.put(Math.abs(ita.getStart() - it.getStart()), ita);
                }
            }

            if (cands.size() > 0) {
                int i = cands.keySet().iterator().next();
                return cands.get(i);
            }
        }

        return null;
    }

    public List<List<Interval>> alignSmoothly(String sFw) {
        List<List<Interval>> kfw = new ArrayList<>();
        List<List<Interval>> krc = new ArrayList<>();

        int ufw = 0, urc = 0;

        for (int i = 0; i <= sFw.length() - kmerSize; i++) {
            String fw = sFw.substring(i, i + kmerSize);
            String rc = SequenceUtils.reverseComplement(fw);

            kfw.add(new ArrayList<>(findKmer(fw)));
            krc.add(flipStrands(new ArrayList<>(findKmer(rc))));

            if (kfw.get(i).size() == 1) { ufw++; }
            if (krc.get(i).size() == 1) { urc++; }
        }

        for (int i = 0; i < kfw.size(); i++) {
            if (kfw.get(i).size() > 1) {
                Interval ita = closestMatchingAlignment(closestUniqueAlignment(kfw, i), kfw.get(i));
                List<Interval> its = new ArrayList<>();

                if (ita != null) { its.add(ita); }

                kfw.set(i, its);
            }

            if (krc.get(i).size() > 1) {
                Interval ita = closestMatchingAlignment(closestUniqueAlignment(krc, i), krc.get(i));
                List<Interval> its = new ArrayList<>();

                if (ita != null) { its.add(ita); }

                krc.set(i, its);
            }
        }


        List<List<Interval>> finalLis = (ufw > urc) ? kfw : krc;

        /*
        List<Interval> finalLi = new ArrayList<Interval>();

        for (List<Interval> li : finalLis) {
            if (li.size() == 1) {
                finalLi.add(li.get(0));
            } else {
                finalLi.add(null);
            }
        }

        return finalLi;
        */

        return finalLis;
    }
}
