package uk.ac.ox.well.cortexjdk.utils.io.graph.cortex;

import com.google.common.base.Joiner;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.utils.LineReader;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;

import java.io.File;
import java.util.*;

public class CortexCollection implements uk.ac.ox.well.cortexjdk.utils.io.graph.DeBruijnGraph {
    private CortexRecord[] nextRecs;

    private List<CortexGraph> graphList = new ArrayList<>();
    private Map<CortexGraph, Pair<List<Integer>, List<Integer>>> graphs = new HashMap<>();
    private Map<CortexGraph, Pair<List<Integer>, List<Integer>>> fgraphs = new HashMap<>();

    private int numColors = 0;
    private int kmerSize = 0;
    private int kmerBits = 0;

    private List<CortexColor> colors = new ArrayList<>();

    public CortexCollection(List<CortexGraph> graphCollection) {
        loadCollection(graphCollection);
    }

    public CortexCollection(CortexGraph... graphCollection) {
        loadCollection(Arrays.asList(graphCollection));
    }

    private void loadCollection(List<CortexGraph> graphCollection) {
        int ac = 0;

        for (CortexGraph g : graphCollection) {
            if (kmerSize == 0) {
                kmerSize = g.getKmerSize();
                kmerBits = g.getKmerBits();
            }

            if (kmerSize != g.getKmerSize()) {
                throw new CortexJDKException("Graph kmer sizes are not equal.  Expected k=" + kmerSize + ", but found k=" + g.getKmerSize() + " in graph " + g.getFile().getAbsolutePath());
            }

            List<Integer> accessColors = new ArrayList<>();
            List<Integer> loadingColors = new ArrayList<>();

            for (int lc = 0; lc < g.getNumColors(); lc++, ac++) {
                accessColors.add(ac);
                loadingColors.add(lc);
            }

            graphs.put(g, new Pair<>(accessColors, loadingColors));
            fgraphs.put(new CortexGraph(g.getFile()), new Pair<>(accessColors, loadingColors));
            graphList.add(g);

            numColors += accessColors.size();

            colors.addAll(g.getColors());
        }
    }

    public CortexGraph getGraph(int color) {
        for (CortexGraph g : graphs.keySet()) {
            Pair<List<Integer>, List<Integer>> p = graphs.get(g);

            if (p.getFirst().contains(color)) {
                return g;
            }
        }

        throw new CortexJDKException("Color doesn't exist in graph.");
    }

    public String getSampleName(int color) {
        CortexGraph g = getGraph(color);

        Pair<List<Integer>, List<Integer>> p = graphs.get(g);

        List<Integer> accessColors = p.getFirst();
        List<Integer> loadingColors = p.getSecond();

        for (int i = 0; i < accessColors.size(); i++) {
            if (accessColors.get(i) == color) {
                CortexColor cc = g.getColor(loadingColors.get(i));

                return cc.getSampleName();
            }
        }

        throw new CortexJDKException("Color " + color + " not found in composite graph");
    }

    public int getNumColors() { return numColors; }

    @Override
    public long getNumRecords() {
        return 0;
    }

    @Override
    public List<CortexColor> getColors() {
        return colors;
    }

    @Override
    public boolean hasColor(int color) {
        return colors.size() > color;
    }

    @Override
    public CortexColor getColor(int color) {
        return colors.get(color);
    }

    @Override
    public int getColorForSampleName(String sampleName) {
        List<Integer> l = getColorsForSampleNames(Collections.singletonList(sampleName));

        return l.size() == 1 ? l.get(0) : -1;
    }

    @Override
    public List<Integer> getColorsForSampleNames(Collection<String> sampleNames) {
        List<Integer> l = new ArrayList<>();

        Set<String> names = new LinkedHashSet<>(sampleNames);

        for (int c = 0; c < colors.size(); c++) {
            if (names.contains(colors.get(c).getSampleName())) {
                l.add(c);
            }
        }

        return l;
    }

    public int getKmerSize() { return kmerSize; }

    public int getKmerBits() { return kmerBits; }

    public CortexHeader getHeader() {
        CortexHeader ch = new CortexHeader();
        ch.setVersion(6);
        ch.setNumColors(getNumColors());
        ch.setKmerSize(getKmerSize());
        ch.setKmerBits(getKmerBits());

        for (int c = 0; c < getNumColors(); c++) {
            ch.addColor(colors.get(c));
        }

        return ch;
    }

    @Override
    public int getVersion() { return 6; }

    public CortexRecord findRecord(byte[] bk) {
        long[] binaryKmer = null;
        int[] coverages = new int[numColors];
        byte[] edges = new byte[numColors];

        for (int assignmentColor = 0; assignmentColor <= numColors; assignmentColor++) {
            for (CortexGraph g : fgraphs.keySet()) {
                CortexRecord cr = g.findRecord(bk);

                if (cr != null) {
                    binaryKmer = cr.getBinaryKmer();

                    Pair<List<Integer>, List<Integer>> p = fgraphs.get(g);
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

    public CortexRecord findRecord(CortexByteKmer bk) {
        return findRecord(bk.getKmer());
    }
    public CortexRecord findRecord(CanonicalKmer ck) {
        return findRecord(ck.getKmerAsBytes());
    }
    public CortexRecord findRecord(String kmer) {
        return findRecord(kmer.getBytes());
    }

    @Override
    public File getFile() {
        throw new UnsupportedOperationException();
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
    public long position() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void position(long i) {
        throw new UnsupportedOperationException();
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

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void close() {
        for (CortexGraph g : graphList) {
            g.close();
        }
    }

    @Override
    public CortexRecord getRecord(long i) {
        throw new UnsupportedOperationException();
    }

    @Override
    public String toString() {
        List<String> filenames = new ArrayList<>();
        for (CortexGraph g : graphList) {
            filenames.add(g.getFile().getAbsolutePath());
        }

        String info = "files:\n  " + Joiner.on("\n  ").join(filenames) + "\n"
                + "----" + "\n"
                + "binary version: " + this.getVersion() + "\n"
                + "kmer size: " + this.getKmerSize() + "\n"
                + "bitfields: " + this.getKmerBits() + "\n"
                + "colors: " + this.getNumColors() + "\n";

        for (int color = 0; color < this.getNumColors(); color++) {
            CortexColor cortexColor = this.getColors().get(color);
            info += "-- Color " + color + " --\n"
                 +  "  sample name: '" + cortexColor.getSampleName() + "'\n"
                 +  "  mean read length: " + cortexColor.getMeanReadLength() + "\n"
                 +  "  total sequence loaded: " + "(not parsed)" + "\n"
                 +  "  sequence error rate: " + "(not parsed)" + "\n"
                 +  "  tip clipping: " + (cortexColor.isTipClippingApplied() ? "yes" : "no") + "\n"
                 +  "  remove_low_coverage_supernodes: " + (cortexColor.isLowCovgSupernodesRemoved() ? "yes" : "no") + "\n"
                 +  "  remove_low_coverage_kmers: " + (cortexColor.isLowCovgKmersRemoved() ? "yes" : "no") + "\n"
                 +  "  cleaned against graph: " + (cortexColor.isCleanedAgainstGraph() ? "yes" : "no") + "\n";
        }

        return info;
    }
}
