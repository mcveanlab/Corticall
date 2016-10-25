package uk.ac.ox.well.indiana.utils.io.cortex.collection;

import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;

import java.io.File;
import java.util.*;

public class CortexCollection2 implements Iterable<CortexRecord>, Iterator<CortexRecord> {
    private Map<Integer, CortexGraph> colorToGraphMap = new TreeMap<>();
    private Map<Integer, Integer> colorToColorMap = new HashMap<>();
    private int maxAssignmentColor = -1;

    private CortexRecord[] nextRecs;

    public CortexCollection2(String collectionFileString) {
        if (collectionFileString.endsWith(".ctx")) {
            CortexGraph cg = new CortexGraph(collectionFileString);

            for (int color = 0; color < cg.getNumColors(); color++) {
                colorToGraphMap.put(color, new CortexGraph(collectionFileString));
                colorToColorMap.put(color, color);
            }

            maxAssignmentColor = cg.getNumColors() - 1;
            cg.close();
        } else {
            LineReader lr = new LineReader(new File(collectionFileString));

            String l;
            while ((l = lr.getNextRecord()) != null) {
                l = l.replaceAll("^\\s+", "");
                l = l.replaceAll("\\s+$", "");

                String[] p = l.split("[\\s:]+");

                int assignmentColor = Integer.valueOf(p[0]);
                String graphFile = p[1];
                int loadingColor = p.length > 2 ? Integer.valueOf(p[2]) : 0;

                CortexGraph cg = new CortexGraph(graphFile);
                colorToGraphMap.put(assignmentColor, cg);
                colorToColorMap.put(assignmentColor, loadingColor);

                if (assignmentColor > maxAssignmentColor) {
                    maxAssignmentColor = assignmentColor;
                }
            }
        }
    }

    public CortexGraph getGraph(int color) {
        return colorToGraphMap.get(color);
    }

    public int getNumColors() { return maxAssignmentColor + 1; }
    public int getKmerSize() { return colorToGraphMap.values().iterator().next().getKmerSize(); }
    public int getKmerBits() { return colorToGraphMap.values().iterator().next().getKmerBits(); }

    public CortexRecord findRecord(String kmer) {
        long[] binaryKmer = null;
        int kmerBits = 0;
        int kmerSize = 0;
        int[] coverages = new int[maxAssignmentColor + 1];
        byte[] edges = new byte[maxAssignmentColor + 1];

        for (int assignmentColor = 0; assignmentColor <= maxAssignmentColor; assignmentColor++) {
            if (colorToGraphMap.containsKey(assignmentColor)) {
                CortexGraph cg = colorToGraphMap.get(assignmentColor);
                int graphColor = colorToColorMap.get(assignmentColor);

                CortexRecord cr = cg.findRecord(kmer);

                if (cr != null) {
                    if (binaryKmer == null) {
                        binaryKmer = cr.getBinaryKmer();
                        kmerBits = cr.getKmerBits();
                        kmerSize = cr.getKmerSize();
                    }

                    coverages[assignmentColor] = cr.getCoverage(graphColor);
                    edges[assignmentColor] = cr.getEdges()[graphColor];
                }
            } else {
                coverages[assignmentColor] = 0;
                edges[assignmentColor] = 0;
            }
        }

        return (binaryKmer == null) ? null : new CortexRecord(binaryKmer, coverages, edges, kmerSize, kmerBits);
    }

    private void moveToBeginningOfRecordsSection() {
        nextRecs = new CortexRecord[maxAssignmentColor + 1];

        for (int color = 0; color <= maxAssignmentColor; color++) {
            nextRecs[color] = colorToGraphMap.get(color) == null ? null : colorToGraphMap.get(color).next();
        }
    }

    @Override
    public Iterator<CortexRecord> iterator() {
        moveToBeginningOfRecordsSection();

        return this;
    }

    @Override
    public boolean hasNext() {
        for (int color = 0; color <= maxAssignmentColor; color++) {
            if (nextRecs[color] != null) {
                return true;
            }
        }

        return false;
    }

    @Override
    public CortexRecord next() {
        int lc = -1;

        for (int color = 0; color <= maxAssignmentColor; color++) {
            if (nextRecs[color] != null) {
                if (lc == -1 || nextRecs[color].getKmerAsString().compareTo(nextRecs[lc].getKmerAsString()) < 0) {
                    lc = color;
                }
            }
        }

        if (lc >= 0) {
            String compKmer = nextRecs[lc].getKmerAsString();

            long[] binaryKmer = null;
            int kmerBits = 0;
            int kmerSize = 0;
            int[] coverages = new int[maxAssignmentColor + 1];
            byte[] edges = new byte[maxAssignmentColor + 1];

            for (int color = 0; color <= maxAssignmentColor; color++) {
                if (nextRecs[color] != null && nextRecs[color].getKmerAsString().equals(compKmer)) {
                    CortexRecord cr = nextRecs[color];

                    if (binaryKmer == null) {
                        binaryKmer = cr.getBinaryKmer();
                        kmerBits = cr.getKmerBits();
                        kmerSize = cr.getKmerSize();
                    }

                    coverages[color] = cr.getCoverage(colorToColorMap.get(color));
                    edges[color] = cr.getEdges()[colorToColorMap.get(color)];

                    CortexRecord ncr = colorToGraphMap.get(color) == null ? null : colorToGraphMap.get(color).next();

                    if (ncr != null && cr.getKmerAsString().compareTo(ncr.getKmerAsString()) >= 0) {
                        throw new IndianaException("Records are not sorted ('" + cr.getKmerAsString() + "' is found before '" + ncr.getKmerAsString() + "' but is lexicographically greater)");
                    }

                    nextRecs[color] = ncr;
                } else {
                    coverages[color] = 0;
                    edges[color] = 0;
                }
            }

            return (binaryKmer == null) ? null : new CortexRecord(binaryKmer, coverages, edges, kmerSize, kmerBits);
        }

        return null;
    }
}
