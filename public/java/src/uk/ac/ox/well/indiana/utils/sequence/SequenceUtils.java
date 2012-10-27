package uk.ac.ox.well.indiana.utils.sequence;

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
            }

            rc[sequence.length - 1 - i] = rcBase;
        }

        return rc;
    }

    public static byte[] getAlphanumericallyLowestOrientation(byte[] sequence) {
        byte[] rc = getReverseComplement(sequence);

        String kmerStr = new String(sequence);
        String rcStr = new String(rc);

        return (kmerStr.compareTo(rcStr) < 0) ? sequence : rc;
    }
}
