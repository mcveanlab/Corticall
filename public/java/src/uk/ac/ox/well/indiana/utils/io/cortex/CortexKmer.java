package uk.ac.ox.well.indiana.utils.io.cortex;

import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.util.Arrays;

public class CortexKmer {
    private byte[] kmer;

    public CortexKmer(String kmer) {
        this.kmer = SequenceUtils.alphanumericallyLowestOrientation(kmer.getBytes());
    }

    public CortexKmer(byte[] kmer) {
        this.kmer = SequenceUtils.alphanumericallyLowestOrientation(kmer);
    }

    public CortexKmer(byte[] kmer, boolean kmerIsAlphanumericallyLowest) {
        this.kmer = kmer;

        if (!kmerIsAlphanumericallyLowest) {
            this.kmer = SequenceUtils.alphanumericallyLowestOrientation(kmer);
        }
    }

    public byte[] getKmerAsBytes() {
        return kmer;
    }

    public String getKmerAsString() {
        return new String(kmer);
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(kmer);
    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof CortexKmer) {
            return hashCode() == obj.hashCode();
        } else if (obj instanceof String) {
            return hashCode() == new CortexKmer((String) obj).hashCode();
        } else if (obj instanceof byte[]) {
            return hashCode() == new CortexKmer((byte[]) obj).hashCode();
        }

        return false;
    }

    @Override
    public String toString() {
        return new String(kmer);
    }
}
