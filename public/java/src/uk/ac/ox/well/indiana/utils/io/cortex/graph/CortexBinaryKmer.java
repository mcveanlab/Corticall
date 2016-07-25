package uk.ac.ox.well.indiana.utils.io.cortex.graph;

import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.util.Arrays;

public class CortexBinaryKmer implements Comparable<CortexBinaryKmer> {
    private long[] binaryKmer;

    public CortexBinaryKmer(long[] binaryKmer) {
        this.binaryKmer = binaryKmer;
    }

    public CortexBinaryKmer(byte[] kmer) {
        this.binaryKmer = CortexUtils.encodeBinaryKmer(SequenceUtils.alphanumericallyLowestOrientation(kmer));
    }

    public long[] getBinaryKmer() {
        return binaryKmer;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        CortexBinaryKmer that = (CortexBinaryKmer) o;

        return Arrays.equals(binaryKmer, that.binaryKmer);

    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(binaryKmer);
    }

    @Override
    public int compareTo(CortexBinaryKmer o) {
        if (!Arrays.equals(binaryKmer, o.getBinaryKmer())) {
            for (int b = 0; b < binaryKmer.length; b++) {
                if (binaryKmer[b] != o.getBinaryKmer()[b]) {
                    return (binaryKmer[b] < o.getBinaryKmer()[b]) ? -1 : 1;
                }
            }
        }

        return 0;
    }
}
