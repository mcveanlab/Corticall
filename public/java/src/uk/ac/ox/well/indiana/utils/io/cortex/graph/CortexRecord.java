package uk.ac.ox.well.indiana.utils.io.cortex.graph;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

public class CortexRecord implements Comparable<CortexRecord> {
    private int kmerSize, kmerBits;
    private long[] kmer;
    private int[] coverages;
    private byte[] edges;

    public CortexRecord(long[] binaryKmer, int[] coverages, byte[] edges, int kmerSize, int kmerBits) {
        this.kmer = new long[binaryKmer.length];
        this.coverages = new int[coverages.length];
        this.edges = new byte[edges.length];

        this.kmer = Arrays.copyOf(binaryKmer, binaryKmer.length);
        this.coverages = Arrays.copyOf(coverages, coverages.length);
        this.edges = Arrays.copyOf(edges, edges.length);

        this.kmerSize = kmerSize;
        this.kmerBits = kmerBits;
    }

    private byte binaryNucleotideToChar(long nucleotide) {
        switch ((int) nucleotide) {
            case 0: return 'A';
            case 1: return 'C';
            case 2: return 'G';
            case 3: return 'T';
            default:
                throw new RuntimeException("Nucleotide '" + nucleotide + "' is not a valid binary nucleotide");
        }
    }

    private byte[] decodeBinaryKmer() {
        byte[] rawKmer = new byte[kmerSize];

        long[] binaryKmer = Arrays.copyOf(kmer, kmer.length);

        for (int i = 0; i < binaryKmer.length; i++) {
            binaryKmer[i] = reverse(binaryKmer[i]);
        }

        for (int i = kmerSize - 1; i >= 0; i--) {
            rawKmer[i] = binaryNucleotideToChar(binaryKmer[kmerBits - 1] & 0x3);

            shiftBinaryKmerByOneBase(binaryKmer, kmerBits);
        }

        return rawKmer;
    }

    private void shiftBinaryKmerByOneBase(long[] binaryKmer, int bitfields) {
        for(int i = bitfields - 1; i > 0; i--) {
            binaryKmer[i] >>>= 2;
            binaryKmer[i] |= (binaryKmer[i-1] << 62); // & 0x3
        }
        binaryKmer[0] >>>= 2;
    }

    private long reverse(long x) {
        ByteBuffer bbuf = ByteBuffer.allocate(8);
        bbuf.order(ByteOrder.BIG_ENDIAN);
        bbuf.putLong(x);
        bbuf.order(ByteOrder.LITTLE_ENDIAN);

        return bbuf.getLong(0);
    }

    public int getKmerSize() { return kmerSize; }
    public int getKmerBits() { return kmerBits; }
    public int getNumColors() { return coverages.length; }

    public long[] getKmer() { return this.kmer; }
    public byte[] getKmerAsBytes() { return decodeBinaryKmer(); }
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
        return Arrays.hashCode(kmer) - Arrays.hashCode(coverages) + Arrays.hashCode(edges);
    }

    public int compareTo(CortexRecord cortexRecord) {
        return getKmerAsString().compareTo(cortexRecord.getKmerAsString());
    }

    public Collection<Byte> getInEdgesAsBytes(int color) {
        Collection<Byte> leftEdges = new ArrayList<Byte>();

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
        Collection<Byte> leftEdges = new ArrayList<Byte>();

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
        Collection<String> leftEdgesAsStrings = new ArrayList<String>();

        for (Byte e : leftEdges) {
            byte[] edge = { e };
            leftEdgesAsStrings.add(new String(edge));
        }

        return leftEdgesAsStrings;
    }

    public Collection<String> getInEdgesComplementAsStrings(int color) {
        Collection<Byte> leftEdges = getInEdgesComplementAsBytes(color);
        Collection<String> leftEdgesAsStrings = new ArrayList<String>();

        for (Byte e : leftEdges) {
            byte[] edge = { e };
            leftEdgesAsStrings.add(new String(edge));
        }

        return leftEdgesAsStrings;
    }

    public Collection<Byte> getOutEdgesAsBytes(int color) {
        Collection<Byte> rightEdges = new ArrayList<Byte>();

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
        Collection<Byte> rightEdges = new ArrayList<Byte>();

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
        Collection<String> rightEdgesAsStrings = new ArrayList<String>();

        for (Byte e : rightEdges) {
            byte[] edge = { e };
            rightEdgesAsStrings.add(new String(edge));
        }

        return rightEdgesAsStrings;
    }

    public Collection<String> getOutEdgesComplementAsStrings(int color) {
        Collection<Byte> rightEdges = getOutEdgesAsBytes(color);
        Collection<String> rightEdgesAsStrings = new ArrayList<String>();

        for (Byte e : rightEdges) {
            byte[] edge = { e };
            rightEdgesAsStrings.add(new String(edge));
        }

        return rightEdgesAsStrings;
    }
}
