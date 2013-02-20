package uk.ac.ox.well.indiana.utils.sequence;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SequenceUtils {
    private SequenceUtils() {}

    public static byte[] getReverseComplement(byte[] sequence) {
        byte[] rc = new byte[sequence.length];

        for (int i = 0; i < sequence.length; i++) {
            byte rcBase = 'N';
            switch (sequence[i]) {
                case 'A': rcBase = 'T'; break;
                case 'a': rcBase = 't'; break;

                case 'C': rcBase = 'G'; break;
                case 'c': rcBase = 'g'; break;

                case 'G': rcBase = 'C'; break;
                case 'g': rcBase = 'c'; break;

                case 'T': rcBase = 'A'; break;
                case 't': rcBase = 'a'; break;

                case 'N': rcBase = 'N'; break;
                case 'n': rcBase = 'n'; break;

                case '.': rcBase = '.'; break;
            }

            rc[sequence.length - 1 - i] = rcBase;
        }

        return rc;
    }

    public static String getReverseComplement(String sequence) {
        return new String(getReverseComplement(sequence.getBytes()));
    }

    public static byte[] getAlphanumericallyLowestOrientation(byte[] sequence) {
        byte[] rc = getReverseComplement(sequence);

        String kmerStr = new String(sequence);
        String rcStr = new String(rc);

        return (kmerStr.compareTo(rcStr) < 0) ? sequence : rc;
    }

    public static String getAlphanumericallyLowestOrientation(String sequence) {
        return new String(getAlphanumericallyLowestOrientation(sequence.getBytes()));
    }

    public static Map<Integer, String> loadSequenceCodesAsAlphanumericallyLowestKmers(List<FastaSequenceFile> fastas, int kmerSize) {
        Map<Integer, String> kmerHash = new HashMap<Integer, String>();

        for (FastaSequenceFile fasta : fastas) {
            kmerHash.putAll(loadSequenceCodesAsAlphanumericallyLowestKmers(fasta, kmerSize));
        }

        return kmerHash;
    }

    public static Map<Integer, String> loadSequenceCodesAsAlphanumericallyLowestKmers(FastaSequenceFile fasta, int kmerSize) {
        Map<Integer, String> kmerHash = new HashMap<Integer, String>();

        ReferenceSequence seq;
        while ((seq = fasta.nextSequence()) != null) {
            kmerHash.putAll(loadSequenceCodesAsAlphanumericallyLowestKmers(seq, kmerSize));
        }

        return kmerHash;
    }

    public static Map<Integer, String> loadSequenceCodesAsAlphanumericallyLowestKmers(ReferenceSequence seq, int kmerSize) {
        Map<Integer, String> kmerHash = new HashMap<Integer, String>();

        byte[] contig = seq.getBases();
        String[] name = seq.getName().split("\\s+");

        for (int i = 0; i < contig.length - kmerSize; i++) {
            String kmer = new String(getAlphanumericallyLowestOrientation(Arrays.copyOfRange(contig, i, i + kmerSize)));

            kmerHash.put(kmer.hashCode(), name[0]);
        }

        return kmerHash;
    }

    public static int numSegregatingSites(String s1, String s2) {
        if (s1.length() != s2.length()) {
            throw new RuntimeException("Cannot compute number of segreating sites between sequences of different lengths (" + s1.length() + " vs " + s2.length() + ")");
        }

        int S = 0;
        for (int i = 0; i < s1.length(); i++) {
            S += (s1.charAt(i) != s2.charAt(i)) ? 1 : 0;
        }

        return S;
    }

}
