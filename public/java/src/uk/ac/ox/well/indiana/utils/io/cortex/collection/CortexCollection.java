package uk.ac.ox.well.indiana.utils.io.cortex.collection;

import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;

import java.io.File;
import java.util.*;

public class CortexCollection {
    private Map<Integer, CortexGraph> colorToGraphMap = new TreeMap<>();
    private Map<Integer, Integer> colorToColorMap = new HashMap<>();
    private int maxAssignmentColor = -1;

    public CortexCollection(File collectionFile) {
        if (collectionFile.getAbsolutePath().endsWith(".ctx")) {
            CortexGraph cg = new CortexGraph(collectionFile);

            for (int color = 0; color < cg.getNumColors(); color++) {
                colorToGraphMap.put(color, new CortexGraph(collectionFile));
                colorToColorMap.put(color, color);
            }

            maxAssignmentColor = cg.getNumColors() - 1;
            cg.close();
        } else {
            LineReader lr = new LineReader(collectionFile);

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
                    binaryKmer = cr.getBinaryKmer();
                    kmerBits = cr.getKmerBits();
                    kmerSize = cr.getKmerSize();

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
}
