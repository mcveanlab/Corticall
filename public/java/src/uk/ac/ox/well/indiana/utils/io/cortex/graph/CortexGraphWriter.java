package uk.ac.ox.well.indiana.utils.io.cortex.graph;

import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;

public class CortexGraphWriter {
    private File cortexFile;
    private FileOutputStream fos;
    private FileChannel channel;

    private final int version = 6;
    private int kmerSize = 0;
    private int kmerBits = 0;
    private int numColors = 0;

    public CortexGraphWriter(File cortexFile) {
        this.cortexFile = cortexFile;
    }

    private void initialize(int kmerSize, int kmerBits, int numColors) {
        this.kmerSize = kmerSize;
        this.kmerBits = kmerBits;
        this.numColors = numColors;

        try {
            fos = new FileOutputStream(cortexFile);
            channel = fos.getChannel();

            ByteBuffer bb = ByteBuffer.allocateDirect(98);
            bb.order(ByteOrder.LITTLE_ENDIAN);
            bb.clear();

            bb.put("CORTEX".getBytes());

            bb.putInt(this.version);
            bb.putInt(this.kmerSize);
            bb.putInt(this.kmerBits);
            bb.putInt(this.numColors);

            for (int color = 0; color < this.numColors; color++) {
                bb.putInt(7504); // mean read length
            }

            for (int color = 0; color < this.numColors; color++) {
                bb.putLong(7504); // total sequence
            }

            for (int color = 0; color < this.numColors; color++) {
                bb.putInt("PF3D7_0115700".length());
                bb.put("PF3D7_0115700".getBytes());
            }

            for (int color = 0; color < this.numColors; color++) {
                byte[] errorRate = new byte[16];
                bb.put(errorRate);
            }

            for (int color = 0; color < this.numColors; color++) {
                bb.put((byte) 0);
                bb.put((byte) 0);
                bb.put((byte) 0);
                bb.put((byte) 0);
                bb.putInt(0);
                bb.putInt(0);

                bb.putInt("undefined".length());
                bb.put("undefined".getBytes());
            }

            bb.put("CORTEX".getBytes());

            bb.flip();

            channel.write(bb);
        } catch (FileNotFoundException e) {
            throw new IndianaException("Unable to open file '" + cortexFile.getAbsolutePath() + "'", e);
        } catch (IOException e) {
            throw new IndianaException("Unable to write header to file '" + cortexFile.getAbsolutePath() + "'", e);
        }
    }

    public void addRecord(CortexRecord record) {
        if (kmerSize == 0) {
            initialize(record.getKmerSize(), record.getKmerBits(), record.getNumColors());
        }

        ByteBuffer bb = ByteBuffer.allocateDirect(8 * this.numColors*4 + this.numColors);
        bb.clear();

        bb.order(ByteOrder.BIG_ENDIAN);

        long[] binaryKmer = record.getKmer();
        for (int i = 0; i < this.kmerBits; i++) {
            System.out.println(record);
            bb.putLong(binaryKmer[i]);
        }

        bb.order(ByteOrder.LITTLE_ENDIAN);

        int[] coverage = record.getCoverages();
        for (int c = 0; c < this.numColors; c++) {
            bb.putInt(coverage[c]);
        }

        byte[] edges = record.getEdges();
        for (int c = 0; c < this.numColors; c++) {
            bb.put(edges[c]);
        }

        bb.flip();

        try {
            channel.write(bb);
        } catch (IOException e) {
            throw new IndianaException("Unable to write record to file '" + cortexFile.getAbsolutePath() + "'", e);
        }
    }

    public void close() {
        try {
            channel.close();
            fos.close();
        } catch (IOException e) {
            throw new IndianaException("Unable to close '" + cortexFile.getAbsolutePath() + "'", e);
        }
    }
}
