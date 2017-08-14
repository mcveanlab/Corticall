package uk.ac.ox.well.cortexjdk.utils.io.cortex.graph;

import org.jetbrains.annotations.NotNull;

import java.util.Arrays;

/**
 * Created by kiran on 09/08/2017.
 */
public class ByteKmer implements Comparable<ByteKmer> {
    private byte[] kmer;

    public ByteKmer(byte[] kmer) { this.kmer = kmer; }

    public ByteKmer(String kmer) { this.kmer = kmer.getBytes(); }

    /*
    public ByteKmer(byte[] s0, byte[] s1) {
        kmer = new byte[s0.length + s1.length];

        System.arraycopy(s0, 0, kmer, 0, s0.length);
        System.arraycopy(s1, 0, kmer, s0.length, s1.length);
    }
    */

    public int length() { return kmer.length; }

    public byte[] getKmer() { return kmer; }

    public void setKmer(byte[] kmer) { this.kmer = kmer; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ByteKmer kmer1 = (ByteKmer) o;

        return Arrays.equals(kmer, kmer1.kmer);

    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(kmer);
    }

    @Override
    public int compareTo(@NotNull ByteKmer o) {
        byte[] okmer = o.getKmer();
        for (int i = 0; i < kmer.length; i++) {
            if (kmer[i] < okmer[i]) { return -1; }
            if (kmer[i] > okmer[i]) { return 1; }
        }

        return 0;
    }

    @Override
    public String toString() {
        return "ByteKmer{" +
                "kmer=" + (new String(kmer)) +
                '}';
    }
}
