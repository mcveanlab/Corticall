package uk.ac.ox.well.indiana.utils.io.cortex.collection;

import com.google.common.base.Joiner;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexColor;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;

import java.io.File;
import java.util.*;

public class CortexCollection implements Iterable<CortexRecord>, Iterator<CortexRecord> {
    private CortexRecord[] nextRecs;

    private List<CortexGraph> graphList = new ArrayList<>();
    private Map<CortexGraph, Pair<List<Integer>, List<Integer>>> graphs = new HashMap<>();

    private int numColors = 0;
    private int kmerSize = 0;
    private int kmerBits = 0;

    private Map<Integer, String> colorList = new HashMap<>();

    public CortexCollection(String collectionFileString) {
        if (!collectionFileString.endsWith(".ctx") && new File(collectionFileString).exists()) {
            List<String> lines = new ArrayList<>();

            LineReader lr = new LineReader(new File(collectionFileString));
            String l;
            while ((l = lr.getNextRecord()) != null) {
                lines.add(l);
            }

            collectionFileString = Joiner.on(";").join(lines);
        }

        if (collectionFileString.contains(";")) {
            String[] pieces = collectionFileString.split(";");

            int ac = 0;

            for (String piece : pieces) {
                String[] parts = piece.split(":");

                String accessMap = "*";
                CortexGraph graph = null;
                String loadingMap = "*";

                for (int i = 0; i < parts.length; i++) {
                    if (new File(parts[i]).exists()) {
                        graph = new CortexGraph(parts[i]);
                    } else {
                        if (graph == null) {
                            accessMap = parts[i];
                        } else {
                            loadingMap = parts[i];
                        }
                    }
                }

                if (graph == null) { throw new IndianaException("No graph specified"); }

                List<Integer> accessColors = new ArrayList<>();
                List<Integer> loadingColors = new ArrayList<>();

                if (loadingMap.equals("*")) {
                    for (int lc = 0; lc < graph.getNumColors(); lc++, ac++) {
                        accessColors.add(ac);
                        loadingColors.add(lc);
                    }
                } else {
                    String[] accessMapList = accessMap.split(",");
                    String[] loadingMapList = loadingMap.split(",");

                    if (accessMap.equals("*")) {
                        for (int i = 0; i < loadingMapList.length; i++, ac++) {
                            accessColors.add(ac);
                            loadingColors.add(Integer.valueOf(loadingMapList[i]));
                        }
                    } else {
                        if (accessMapList.length != loadingMapList.length) {
                            throw new IndianaException("Color access/loading lists are unequal length");
                        }

                        for (int i = 0; i < accessMapList.length; i++) {
                            ac = Integer.valueOf(accessMapList[i]);
                            int lc = Integer.valueOf(loadingMapList[i]);

                            accessColors.add(ac);
                            loadingColors.add(lc);
                        }
                    }
                }

                graphs.put(graph, new Pair<>(accessColors, loadingColors));
                graphList.add(graph);

                numColors += accessColors.size();
            }
        } else if (collectionFileString.endsWith(".ctx")) {
            CortexGraph graph = new CortexGraph(collectionFileString);

            List<Integer> accessColors = new ArrayList<>();
            List<Integer> loadingColors = new ArrayList<>();

            for (int c = 0; c < graph.getNumColors(); c++) {
                accessColors.add(c);
                loadingColors.add(c);
            }

            graphs.put(graph, new Pair<>(accessColors, loadingColors));
            graphList.add(graph);

            numColors = graph.getNumColors();
        }

        for (CortexGraph g : graphs.keySet()) {
            if (kmerSize == 0) {
                kmerSize = g.getKmerSize();
            } else {
                if (g.getKmerSize() != kmerSize) {
                    throw new IndianaException("Kmer size between graphs don't match.");
                }
            }

            if (kmerBits == 0) {
                kmerBits = g.getKmerBits();
            } else {
                if (g.getKmerBits() != kmerBits) {
                    throw new IndianaException("Kmer bits between graphs don't match.");
                }
            }
        }
    }

    public CortexGraph getGraph(int color) {
        for (CortexGraph g : graphs.keySet()) {
            Pair<List<Integer>, List<Integer>> p = graphs.get(g);

            if (p.getFirst().contains(color)) {
                return g;
            }
        }

        throw new IndianaException("Color doesn't exist in graph.");
    }

    public String getSampleName(int color) {
        CortexGraph g = getGraph(color);

        Pair<List<Integer>, List<Integer>> p = graphs.get(g);

        List<Integer> accessColors = p.getFirst();
        List<Integer> loadingColors = p.getSecond();

        for (int i = 0; i < accessColors.size(); i++) {
            if (i == color) {
                CortexColor cc = g.getColor(loadingColors.get(i));

                return cc.getSampleName();
            }
        }

        throw new IndianaException("Color " + color + " not found in composite graph");
    }

    public int getNumColors() { return numColors; }

    public int getKmerSize() { return kmerSize; }

    public int getKmerBits() { return kmerBits; }

    public CortexRecord findRecord(String kmer) {
        long[] binaryKmer = null;
        int[] coverages = new int[numColors];
        byte[] edges = new byte[numColors];

        for (int assignmentColor = 0; assignmentColor <= numColors; assignmentColor++) {
            for (CortexGraph g : graphs.keySet()) {
                CortexRecord cr = g.findRecord(kmer);

                if (cr != null) {
                    binaryKmer = cr.getBinaryKmer();

                    Pair<List<Integer>, List<Integer>> p = graphs.get(g);
                    List<Integer> accessColors = p.getFirst();
                    List<Integer> loadingColors = p.getSecond();

                    for (int c = 0; c < accessColors.size(); c++) {
                        int ac = accessColors.get(c);
                        int lc = loadingColors.get(c);

                        coverages[ac] = cr.getCoverage(lc);
                        edges[ac] = cr.getEdges()[lc];
                    }
                }
            }
        }

        return (binaryKmer == null) ? null : new CortexRecord(binaryKmer, coverages, edges, kmerSize, kmerBits);
    }

    private void moveToBeginningOfRecordsSection() {
        nextRecs = new CortexRecord[graphs.keySet().size() + 1];

        int i = 0;
        for (CortexGraph g : graphList) {
            nextRecs[i] = g.next();
            i++;
        }
    }

    @Override
    public Iterator<CortexRecord> iterator() {
        moveToBeginningOfRecordsSection();

        return this;
    }

    @Override
    public boolean hasNext() {
        for (int i = 0; i < nextRecs.length; i++) {
            if (nextRecs[i] != null) {
                return true;
            }
        }

        return false;
    }

    @Override
    public CortexRecord next() {
        int lowc = -1;

        for (int i = 0; i < nextRecs.length; i++) {
            if (nextRecs[i] != null) {
                if (lowc == -1 || nextRecs[i].getKmerAsString().compareTo(nextRecs[lowc].getKmerAsString()) < 0) {
                    lowc = i;
                }
            }
        }

        if (lowc >= 0) {
            String compKmer = nextRecs[lowc].getKmerAsString();

            long[] binaryKmer = null;
            int[] coverages = new int[numColors];
            byte[] edges = new byte[numColors];

            for (int i = 0; i < nextRecs.length; i++) {
                if (nextRecs[i] != null && nextRecs[i].getKmerAsString().equals(compKmer)) {
                    CortexRecord cr = nextRecs[i];

                    if (cr != null) {
                        binaryKmer = cr.getBinaryKmer();

                        CortexGraph g = graphList.get(i);

                        Pair<List<Integer>, List<Integer>> p = graphs.get(g);
                        List<Integer> accessColors = p.getFirst();
                        List<Integer> loadingColors = p.getSecond();

                        for (int c = 0; c < accessColors.size(); c++) {
                            int ac = accessColors.get(c);
                            int lc = loadingColors.get(c);

                            coverages[ac] = cr.getCoverage(lc);
                            edges[ac] = cr.getEdges()[lc];
                        }

                        nextRecs[i] = graphList.get(i).next();
                    }
                }
            }

            return (binaryKmer == null) ? null : new CortexRecord(binaryKmer, coverages, edges, kmerSize, kmerBits);
        }

        return null;
    }
}
