package uk.ac.ox.well.indiana.utils.io.cortex;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Arrays;

public class CortexRecord implements Comparable<CortexRecord> {
    private byte[] rawKmer;
    private int[] coverages;
    private byte[] edges;

    public CortexRecord(long[] binaryKmer, int[] coverages, byte[] edges, int kmerSize, int kmerBits) {
        this.coverages = coverages;
        this.edges = edges;
        this.rawKmer = decodeBinaryKmer(binaryKmer, kmerSize, kmerBits);
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

    private byte[] decodeBinaryKmer(long[] binaryKmer, int kmerSize, int kmerBits) {
        byte[] rawKmer = new byte[kmerSize];

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

    public byte[] getKmer() {
        return rawKmer;
    }

    public String getKmerString() {
        return new String(getKmer());
    }

    public String[] getEdges() {
        int numColors = edges.length;
        byte[] str = {'a', 'c', 'g', 't', 'A', 'C', 'G', 'T'};

        String[] edgesString = new String[numColors];
        for (int color = 0; color < numColors; color++) {
            byte edge = edges[color];

            int left = (edge >> 4);
            int right = (edge & 0xf);

            byte[] edgeStr = new byte[8];

            for (int i = 0; i < 4; i++) {
                int leftEdge = (left & (0x1 << (3-i)));
                edgeStr[i] = (byte) ((leftEdge != 0) ? str[i] : '.');

                int rightEdge = (right & (0x1 << i));
                edgeStr[i+4] = (byte) ((rightEdge != 0) ? str[i+4] : '.');
            }

            edgesString[color] = new String(edgeStr);
        }

        return edgesString;
    }

    public String getEdges(int color) {
        return getEdges()[color];
    }

    public int[] getCoverages() {
        return coverages;
    }

    public int getCoverage(int color) {
        return getCoverages()[color];
    }

    public String toString() {
        String info = getKmerString();

        for (int coverage : getCoverages()) {
            info += " " + coverage;
        }

        for (String edge : getEdges()) {
            info += " " + edge;
        }

        return info;
    }

    public int hashCode() {
        return Arrays.hashCode(rawKmer) - Arrays.hashCode(coverages) + Arrays.hashCode(edges);
    }

    public int compareTo(CortexRecord cortexRecord) {
        return getKmerString().compareTo(cortexRecord.getKmerString());
    }
}
