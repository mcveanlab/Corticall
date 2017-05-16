package uk.ac.ox.well.indiana.utils.io.cortex.graph;

import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
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

        constructRecord(sk, coverageList, inEdgesList, outEdgesList);
    }

    public CortexRecord(String recordString) {
        String[] fields = recordString.split("\\s+");

        String sk = fields[0];
        int numColors = (fields.length - 1) / 2;

        List<Integer> coverageList = new ArrayList<>();
        List<Set<String>> inEdgesList = new ArrayList<>();
        List<Set<String>> outEdgesList = new ArrayList<>();

        for (int c = 0; c < numColors; c++) {
            Set<String> inEdges = new HashSet<>();
            Set<String> outEdges = new HashSet<>();

            String cov = fields[1 + c];
            coverageList.add(Integer.valueOf(cov));

            String edgesString = (fields[1 + numColors + c]).toUpperCase();

            for (int i = 0; i < 4; i++) {
                if (edgesString.charAt(i) != '.') {
                    inEdges.add(String.valueOf(edgesString.charAt(i)));
                }
            }

            for (int i = 4; i < 8; i++) {
                if (edgesString.charAt(i) != '.') {
                    outEdges.add(String.valueOf(edgesString.charAt(i)));
                }
            }

            inEdgesList.add(inEdges);
            outEdgesList.add(outEdges);
        }

        constructRecord(sk, coverageList, inEdgesList, outEdgesList);
    }

    private void constructRecord(String sk, List<Integer> coverageList, List<Set<String>> inEdgesList, List<Set<String>> outEdgesList) {
        CortexKmer ck = new CortexKmer(sk);
        this.kmerSize = ck.length();
        this.kmerBits = getKmerBits(kmerSize);

        int colors = coverageList.size();
        this.binaryKmer = new CortexBinaryKmer(sk.getBytes()).getBinaryKmer();
        this.coverages = new int[colors];
        this.edges = new byte[colors];

        for (int c = 0; c < colors; c++) {
            coverages[c] = coverageList.get(c);

            Set<String> inEdges = inEdgesList.get(c);
            Set<String> outEdges = outEdgesList.get(c);

            byte edge = encodeBinaryEdges(inEdges, outEdges, ck.isFlipped());

            edges[c] = edge;
        }
    }

    public int getKmerSize() { return kmerSize; }
    public int getKmerBits() { return kmerBits; }
    public int getNumColors() { return coverages.length; }

    public long[] getBinaryKmer() { return this.binaryKmer; }
    public byte[] getKmerAsBytes() { return decodeBinaryKmer(binaryKmer, kmerSize, kmerBits); }
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

    public int[] getCoverages() { return coverages; }

    public int getCoverage(int color) { return getCoverages()[color]; }

    public byte[] getEdgesAsBytes(int color) { return getEdgesAsBytes()[color]; }

    public String getEdgesAsString(int color) {
        return getEdgeAsStrings()[color];
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

    public Collection<Byte> getInEdgesAsBytes(int color, boolean complement) {
        Collection<Byte> leftEdges = new ArrayList<>();

        byte[] str = {'A', 'C', 'G', 'T'};

        if (complement) {
            str = new byte[] {'T', 'G', 'C', 'A'};
        }

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

    public Collection<String> getInEdgesAsStrings(int color, boolean complement) {
        Collection<Byte> leftEdges = getInEdgesAsBytes(color, complement);
        Collection<String> leftEdgesAsStrings = new ArrayList<>();

        for (Byte e : leftEdges) {
            byte[] edge = { e };
            leftEdgesAsStrings.add(new String(edge));
        }

        return leftEdgesAsStrings;
    }

    public Collection<Byte> getOutEdgesAsBytes(int color, boolean complement) {
        Collection<Byte> rightEdges = new ArrayList<>();

        byte[] str = {'A', 'C', 'G', 'T'};

        if (complement) {
            str = new byte[] {'T', 'G', 'C', 'A'};
        }

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

    public Collection<String> getOutEdgesAsStrings(int color, boolean complement) {
        Collection<Byte> rightEdges = getOutEdgesAsBytes(color, complement);
        Collection<String> rightEdgesAsStrings = new ArrayList<>();

        for (Byte e : rightEdges) {
            byte[] edge = { e };
            rightEdgesAsStrings.add(new String(edge));
        }

        return rightEdgesAsStrings;
    }

    public int getInDegree(int color) { return getInEdgesAsBytes(color, false).size(); }

    public int getOutDegree(int color) { return getOutEdgesAsBytes(color, false).size(); }

    public static byte[] decodeBinaryKmer(long[] kmer, int kmerSize, int kmerBits) {
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

    public static int getKmerBits(int kmerSize) {
        return (int) Math.ceil(((float) kmerSize)/32.0);
    }

    public static long[] encodeBinaryKmer(byte[] kmer) {
        int numBits = getKmerBits(kmer.length);
        long[] binaryKmer = new long[numBits];

        for (int b = 0; b < numBits; b++) {
            for (int i = kmer.length - 32*(b+1); i < kmer.length - 32*b; i++) {
                if (i >= 0) {
                    long nuc = charToBinaryNucleotide(kmer[i]);

                    binaryKmer[numBits - b - 1] |= nuc;
                }

                if (i < kmer.length - 32*b - 1) {
                    binaryKmer[numBits - b - 1] <<= 2;
                }
            }

            binaryKmer[numBits - b - 1] = reverse(binaryKmer[numBits - b - 1]);
        }

        return binaryKmer;
    }

    private static byte binaryNucleotideToChar(long nucleotide) {
        switch ((int) nucleotide) {
            case 0: return 'A';
            case 1: return 'C';
            case 2: return 'G';
            case 3: return 'T';
            default:
                throw new RuntimeException("Nucleotide '" + nucleotide + "' is not a valid binary nucleotide");
        }
    }

    private static long charToBinaryNucleotide(byte b) {
        switch (b) {
            case 'A' : return 0;
            case 'C' : return 1;
            case 'G' : return 2;
            case 'T' : return 3;
            case 'a' : return 0;
            case 'c' : return 1;
            case 'g' : return 2;
            case 't' : return 3;
            default:
                throw new RuntimeException("Nucleotide '" + b + "' is not a valid character nucleotide");
        }
    }

    private static void shiftBinaryKmerByOneBase(long[] binaryKmer, int bitfields) {
        for(int i = bitfields - 1; i > 0; i--) {
            binaryKmer[i] >>>= 2;
            binaryKmer[i] |= (binaryKmer[i-1] << 62); // & 0x3
        }
        binaryKmer[0] >>>= 2;
    }

    private static long reverse(long x) {
        ByteBuffer bbuf = ByteBuffer.allocate(8);
        bbuf.order(ByteOrder.BIG_ENDIAN);
        bbuf.putLong(x);
        bbuf.order(ByteOrder.LITTLE_ENDIAN);

        return bbuf.getLong(0);
    }

    private static byte encodeBinaryEdges(Set<String> inEdges, Set<String> outEdges, boolean reverseComplement) {
        byte edge = 0;

        String[] alphabetFwd = { "A", "C", "G", "T" };
        String[] alphabetRev = { "T", "G", "C", "A" };

        if (reverseComplement) {
            String[] a = alphabetFwd;
            alphabetFwd = alphabetRev;
            alphabetRev = a;
        }

        for (int i = 0; i < 4; i++) {
            if (inEdges.contains(alphabetFwd[i])) {
                edge |= 1;
            }
            edge <<= 1;
        }

        for (int i = 0; i < 4; i++) {
            if (outEdges.contains(alphabetRev[i])) {
                edge |= 1;
            }
            if (i != 3) {
                edge <<= 1;
            }
        }

        return edge;
    }
}
