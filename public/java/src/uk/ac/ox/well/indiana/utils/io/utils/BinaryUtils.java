package uk.ac.ox.well.indiana.utils.io.utils;

public class BinaryUtils {
    private BinaryUtils() {}

    public static int toUnsignedInt(byte[] b) {
        long l = 0;
        l |= b[3] & 0xFF;
        l <<= 8;
        l |= b[2] & 0xFF;
        l <<= 8;
        l |= b[1] & 0xFF;
        l <<= 8;
        l |= b[0] & 0xFF;

        return (int) l;
    }
}
