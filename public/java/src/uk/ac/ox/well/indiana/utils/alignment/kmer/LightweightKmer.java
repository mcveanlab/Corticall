package uk.ac.ox.well.indiana.utils.alignment.kmer;

import java.util.Arrays;

public class LightweightKmer {
    private byte[] bk;

    public LightweightKmer(String sk) {
        bk = sk.getBytes();
    }

    public LightweightKmer(String s, int start, int stop) {
        bk = s.substring(start, stop).getBytes();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        LightweightKmer that = (LightweightKmer) o;

        return Arrays.equals(bk, that.bk);
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(bk);
    }

    @Override
    public String toString() {
        return new String(bk);
    }
}
