package uk.ac.ox.well.indiana.utils.io.cortex;

public class CortexRecord {
    private long[] binaryKmer;
    private int[] coverages;
    private byte[] edges;

    private int kmerSize;
    private int kmerBits;

    public CortexRecord(long[] binaryKmer, int[] coverages, byte[] edges, int kmerSize, int kmerBits) {
        this.binaryKmer = binaryKmer;
        this.coverages = coverages;
        this.edges = edges;
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

    private void shiftBinaryKmerByOneBase(long[] binaryKmer, int bitfields) {
        int i;
        for(i = bitfields - 1; i > 0; i--) {
            binaryKmer[i] >>= 2;
            binaryKmer[i] |= (binaryKmer[i-1] << 62); // & 0x3
        }
        binaryKmer[0] >>= 2;
    }

    public String getKmer() {
        byte[] kmer = new byte[kmerSize];
        for (int i = kmerSize - 1; i >= 0; i--) {
            kmer[i] = binaryNucleotideToChar(binaryKmer[kmerBits - 1] & 0x3);
            shiftBinaryKmerByOneBase(binaryKmer, kmerBits);
        }
        String kmerString = new String(kmer);

        return kmerString;
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

    public int[] getCoverages() {
        return coverages;
    }

    public String toString() {
        String info = getKmer();

        int[] coverages = getCoverages();
        for (int color = 0; color < coverages.length; color++) {
            info += " " + coverages[color];
        }

        String[] edges = getEdges();
        for (int color = 0; color < edges.length; color++) {
            info += " " + edges[color];
        }

        return info;
    }
}
