package uk.ac.ox.well.indiana.utils.io.cortex.graph;

import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.util.*;

public class CortexRecord implements Comparable<CortexRecord> {
    private int kmerSize, kmerBits;
    private long[] binaryKmer;
    private int[] coverages;
    private byte[] edges;

    public CortexRecord(long[] binaryKmer, int[] coverages, byte[] edges, int kmerSize, int kmerBits) {
        this.binaryKmer = new long[binaryKmer.length];
        this.coverages = new int[coverages.length];
        this.edges = new byte[edges.length];

        this.binaryKmer = Arrays.copyOf(binaryKmer, binaryKmer.length);
        this.coverages = Arrays.copyOf(coverages, coverages.length);
        this.edges = Arrays.copyOf(edges, edges.length);

        this.kmerSize = kmerSize;
        this.kmerBits = kmerBits;
    }

    public CortexRecord(String sk, List<Integer> coverageList, List<Set<String>> inEdgesList, List<Set<String>> outEdgesList) {
        if (coverageList.size() != inEdgesList.size() && coverageList.size() != outEdgesList.size()) {
            throw new IndianaException("Coverage, in-edge, and out-edge lists must be equal length.");
        }

        CortexKmer ck = new CortexKmer(sk);
        this.kmerSize = ck.length();
        this.kmerBits = CortexUtils.getKmerBits(kmerSize);

        int colors = coverageList.size();
        this.binaryKmer = new CortexBinaryKmer(sk.getBytes()).getBinaryKmer();
        this.coverages = new int[colors];
        this.edges = new byte[colors];

        for (int c = 0; c < colors; c++) {
            coverages[c] = coverageList.get(c);

            Set<String> inEdges = inEdgesList.get(c);
            Set<String> outEdges = outEdgesList.get(c);

            byte edge = CortexUtils.encodeBinaryEdges(inEdges, outEdges, ck.isFlipped());

            edges[c] = edge;
        }
    }

    public int getKmerSize() { return kmerSize; }
    public int getKmerBits() { return kmerBits; }
    public int getNumColors() { return coverages.length; }

    public long[] getBinaryKmer() { return this.binaryKmer; }
    public byte[] getKmerAsBytes() { return CortexUtils.decodeBinaryKmer(binaryKmer, kmerSize, kmerBits); }
    public CortexBinaryKmer getCortexBinaryKmer() { return new CortexBinaryKmer(this.binaryKmer); }
    public CortexKmer getCortexKmer() { return new CortexKmer(getKmerAsBytes(), true); }
    public String getKmerAsString() { return getCortexKmer().getKmerAsString(); }

    public byte[] getEdges() { return edges; }

    public byte[][] getEdgesAsBytes() {
        int numColors = edges.length;
        byte[] str = {'a', 'c', 'g', 't', 'A', 'C', 'G', 'T'};

        byte[][] edgesTable = new byte[numColors][8];
        for (int color = 0; color < numColors; color++) {
            byte edge = edges[color];

            int left = (edge >> 4);
            int right = (edge & 0xf);

            for (int i = 0; i < 4; i++) {
                int leftEdge = (left & (0x1 << (3-i)));
                edgesTable[color][i] = (byte) ((leftEdge != 0) ? str[i] : '.');

                int rightEdge = (right & (0x1 << i));
                edgesTable[color][i+4] = (byte) ((rightEdge != 0) ? str[i+4] : '.');
            }
        }

        return edgesTable;
    }

    public String[] getEdgeAsStrings() {
        byte[][] edgesTable = getEdgesAsBytes();
        String[] edgesStrings = new String[edges.length];

        for (int i = 0; i < edges.length; i++) {
            edgesStrings[i] = new String(edgesTable[i]);
        }

        return edgesStrings;
    }

    public byte[] getEdgesAsBytes(int color) {
        return getEdgesAsBytes()[color];
    }

    public String getEdgesAsString(int color) {
        return getEdgeAsStrings()[color];
    }

    public int[] getCoverages() {
        return coverages;
    }

    public int getCoverage(int color) {
        return getCoverages()[color];
    }

    public String toString() {
        String info = getKmerAsString();

        for (int coverage : getCoverages()) {
            info += " " + coverage;
        }

        for (String edge : getEdgeAsStrings()) {
            info += " " + edge;
        }

        return info;
    }

    public int hashCode() {
        return Arrays.hashCode(binaryKmer) - Arrays.hashCode(coverages) + Arrays.hashCode(edges);
    }

    public boolean equals(Object o) {
        if (o instanceof CortexRecord) {
            CortexRecord o1 = ((CortexRecord) o);

            return (Arrays.equals(binaryKmer, o1.binaryKmer) && Arrays.equals(coverages, o1.coverages) && Arrays.equals(edges, o1.edges));
        }

        return false;
    }

    public int compareTo(CortexRecord cortexRecord) {
        return getKmerAsString().compareTo(cortexRecord.getKmerAsString());
    }

    public Collection<Byte> getInEdgesAsBytes(int color) {
        Collection<Byte> leftEdges = new ArrayList<>();

        byte[] str = {'A', 'C', 'G', 'T'};

        byte edge = edges[color];

        int left = (edge >> 4);

        for (int i = 0; i < 4; i++) {
            int leftEdge = (left & (0x1 << (3-i)));
            if (leftEdge != 0) {
                byte edgeByte = str[i];

                leftEdges.add(edgeByte);
            }
        }

        return leftEdges;
    }

    public Collection<Byte> getInEdgesComplementAsBytes(int color) {
        Collection<Byte> leftEdges = new ArrayList<>();

        byte[] str = {'T', 'G', 'C', 'A'};

        byte edge = edges[color];

        int left = (edge >> 4);

        for (int i = 0; i < 4; i++) {
            int leftEdge = (left & (0x1 << (3-i)));
            if (leftEdge != 0) {
                byte edgeByte = str[i];

                leftEdges.add(edgeByte);
            }
        }

        return leftEdges;
    }

    public Collection<String> getInEdgesAsStrings(int color) {
        Collection<Byte> leftEdges = getInEdgesAsBytes(color);
        Collection<String> leftEdgesAsStrings = new ArrayList<>();

        for (Byte e : leftEdges) {
            byte[] edge = { e };
            leftEdgesAsStrings.add(new String(edge));
        }

        return leftEdgesAsStrings;
    }

    public Collection<String> getInEdgesComplementAsStrings(int color) {
        Collection<Byte> leftEdges = getInEdgesComplementAsBytes(color);
        Collection<String> leftEdgesAsStrings = new ArrayList<>();

        for (Byte e : leftEdges) {
            byte[] edge = { e };
            leftEdgesAsStrings.add(new String(edge));
        }

        return leftEdgesAsStrings;
    }

    public Collection<Byte> getOutEdgesAsBytes(int color) {
        Collection<Byte> rightEdges = new ArrayList<>();

        byte[] str = {'A', 'C', 'G', 'T'};

        byte edge = edges[color];

        int right = (edge & 0xf);

        for (int i = 0; i < 4; i++) {
            int rightEdge = (right & (0x1 << i));

            if (rightEdge != 0) {
                rightEdges.add(str[i]);
            }
        }

        return rightEdges;
    }

    public Collection<Byte> getOutEdgesComplementAsBytes(int color) {
        Collection<Byte> rightEdges = new ArrayList<>();

        byte[] str = {'T', 'G', 'C', 'A'};

        byte edge = edges[color];

        int right = (edge & 0xf);

        for (int i = 0; i < 4; i++) {
            int rightEdge = (right & (0x1 << i));

            if (rightEdge != 0) {
                rightEdges.add(str[i]);
            }
        }

        return rightEdges;
    }

    public Collection<String> getOutEdgesAsStrings(int color) {
        Collection<Byte> rightEdges = getOutEdgesAsBytes(color);
        Collection<String> rightEdgesAsStrings = new ArrayList<>();

        for (Byte e : rightEdges) {
            byte[] edge = { e };
            rightEdgesAsStrings.add(new String(edge));
        }

        return rightEdgesAsStrings;
    }

    public Collection<String> getOutEdgesComplementAsStrings(int color) {
        Collection<Byte> rightEdges = getOutEdgesComplementAsBytes(color);
        Collection<String> rightEdgesAsStrings = new ArrayList<>();

        for (Byte e : rightEdges) {
            byte[] edge = { e };
            rightEdgesAsStrings.add(new String(edge));
        }

        return rightEdgesAsStrings;
    }

    public int getInDegree(int color) { return getInEdgesAsBytes(color).size(); }

    public int getOutDegree(int color) { return getOutEdgesAsBytes(color).size(); }
}
