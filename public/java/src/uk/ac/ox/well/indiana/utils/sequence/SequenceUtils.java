package uk.ac.ox.well.indiana.utils.sequence;

public class SequenceUtils {
    private SequenceUtils() {}

    public static byte[] getReverseComplement(byte[] kmer) {
        byte[] rc = new byte[kmer.length];

        for (int i = 0; i < kmer.length; i++) {
            byte rcBase = 'N';
            switch (kmer[i]) {
                case 'A':
                    rcBase = 'T'; break;
                case 'C':
                    rcBase = 'G'; break;
                case 'G':
                    rcBase = 'C'; break;
                case 'T':
                    rcBase = 'A'; break;
            }

            rc[kmer.length - 1 - i] = rcBase;
        }

        return rc;
    }

    public static byte[] getCortexCompatibleOrientation(byte[] kmer) {
        byte[] rc = getReverseComplement(kmer);

        String kmerStr = new String(kmer);
        String rcStr = new String(rc);

        return (kmerStr.compareTo(rcStr) < 0) ? kmer : rc;
    }
}
