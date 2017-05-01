package uk.ac.ox.well.indiana.utils.io.cortex.graph;

import htsjdk.samtools.util.StringUtil;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class CortexIndex {
    private List<Pair<String, Pair<Long, Long>>> idx = new ArrayList<>();

    public CortexIndex(String idxpath, CortexGraph graph) {
        File idxfile = idxpath != null ? new File(idxpath) : null;

        if (idxfile != null && idxfile.exists()) {
            LineReader lr = new LineReader(idxfile);

            // There's an off-by-one bug in McCortex's index creation, which accumulates per record in the index.
            // This means that the start record in the first line of the index is correct, the second is one higher
            // than it should be, the third is two off, the fourth is three off, etc..  We compensate for that here,
            // but really it needs to be fixed in McCortex, at which point, all of this will break again.
            int erroraccumulator = 0;

            String line;
            while ((line = lr.getNextRecord()) != null) {
                if (!line.startsWith("#")) {
                    String[] pieces = line.split("\t");

                    String kmer = pieces[2];
                    long start = Long.valueOf(pieces[3]) - erroraccumulator;
                    long stop = Long.valueOf(pieces[4]) - 2 - erroraccumulator; // this ensures we get a closed interval

                    Pair<String, Pair<Long, Long>> p = new Pair<>(kmer, new Pair<>(start, stop));

                    idx.add(p);

                    erroraccumulator++;
                }
            }

            validateIndex(graph);
        } else {
            Pair<String, Pair<Long, Long>> p;

            if (graph.getNumRecords() > 0) {
                String kmer = graph.getRecord(0).getKmerAsString();

                p = new Pair<>(kmer, new Pair<>(0L, graph.getNumRecords() - 1));
            } else {
                // For an empty graph
                p = new Pair<>(StringUtil.repeatCharNTimes('A', graph.getKmerSize()), new Pair<>(0L, graph.getNumRecords() - 1));
            }

            idx.add(p);
        }
    }

    private void validateIndex(CortexGraph g) {
        for (Pair<String, Pair<Long, Long>> p : idx) {
            String sk = p.getFirst();
            long offset = p.getSecond().getFirst();

            CortexRecord cr = g.getRecord(offset);

            if (!cr.getKmerAsString().equals(sk)) {
                throw new IndianaException("Index conflicts with graph: in file '" + g.getCortexFile().getAbsolutePath() + "', record " + offset + " expected to be '" + sk + "' but instead was '" + cr.getKmerAsString() + "'.");
            }
        }
    }

    public List<Pair<String, Pair<Long, Long>>> getIndex() {
        return idx;
    }

    public Pair<Long, Long> getBounds(String kmer) {
        if (idx.size() == 1) {
            return new Pair<>(idx.get(0).getSecond());
        }

        int startIndex = 0;
        int stopIndex = idx.size() - 1;
        int midIndex = startIndex + (stopIndex - startIndex) / 2;

        while (startIndex != midIndex && midIndex != stopIndex) {
            String startKmer = idx.get(startIndex).getFirst();
            String midKmer = idx.get(midIndex).getFirst();
            String stopKmer = idx.get(stopIndex).getFirst();

            if (kmer.compareTo(startKmer) < 0) { return null; }
            else if (startKmer.equals(kmer)) { return idx.get(startIndex).getSecond(); }
            else if (midKmer.equals(kmer)) { return idx.get(midIndex).getSecond(); }
            else if (stopKmer.equals(kmer)) { return idx.get(stopIndex).getSecond(); }
            else if (startKmer.compareTo(kmer) < 0 && kmer.compareTo(midKmer) < 0) {
                if (midIndex - startIndex == 1) {
                    return idx.get(startIndex).getSecond();
                } else {
                    stopIndex = midIndex;
                    midIndex = startIndex + ((stopIndex - startIndex) / 2);
                }
            } else if (midKmer.compareTo(kmer) < 0 && kmer.compareTo(stopKmer) < 0) {
                if (stopIndex - midIndex == 1) {
                    return idx.get(midIndex).getSecond();
                } else {
                    startIndex = midIndex;
                    midIndex = startIndex + ((stopIndex - startIndex) / 2);
                }
            } else if (stopKmer.compareTo(kmer) < 0) {
                return idx.get(stopIndex).getSecond();
            }
        }

        return null;
    }
}
