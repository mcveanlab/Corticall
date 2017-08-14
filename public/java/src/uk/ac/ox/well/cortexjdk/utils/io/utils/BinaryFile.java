package uk.ac.ox.well.cortexjdk.utils.io.utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.math.BigInteger;

public class BinaryFile extends RandomAccessFile {
    public BinaryFile(String s, String s1) throws FileNotFoundException {
        super(s, s1);
    }

    public BinaryFile(File file, String s) throws FileNotFoundException {
        super(file, s);
    }

    public int readUnsignedInt() throws IOException {
        byte[] b = new byte[4];
        this.read(b);

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

    public long readUnsignedLong() throws IOException {
        BigInteger b = BigInteger.valueOf(readLong());

        return b.longValue();
    }

    public double readUnsignedDouble() throws IOException {
        long l = this.readLong();
        return Double.longBitsToDouble(Long.reverseBytes(l));
    }
}
