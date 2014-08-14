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

    public CortexGraphWriter(String cortexFilePath) {
        this.cortexFile = new File(cortexFilePath);
    }

    private void initialize(CortexRecord record) {
        this.kmerSize = record.getKmerSize();
        this.kmerBits = record.getKmerBits();
        this.numColors = record.getNumColors();

        try {
            fos = new FileOutputStream(cortexFile);
            channel = fos.getChannel();

            int stringLengths = 0;
            for (CortexColor c : record.getParentGraph().getColors()) {
                stringLengths += c.getSampleName().length() + c.getCleanedAgainstGraphName().length();
            }

            ByteBuffer bb = ByteBuffer.allocateDirect(2 * (44 + numColors*32 + stringLengths));
            bb.order(ByteOrder.LITTLE_ENDIAN);
            bb.clear();

            bb.put("CORTEX".getBytes());

            bb.putInt(version);
            bb.putInt(kmerSize);
            bb.putInt(kmerBits);
            bb.putInt(numColors);

            for (int color = 0; color < numColors; color++) {
                int meanReadLength = record.getParentGraph().getColor(color).getMeanReadLength();
                bb.putInt(meanReadLength);
            }

            for (int color = 0; color < numColors; color++) {
                long totalSequence = record.getParentGraph().getColor(color).getTotalSequence();

                bb.order(ByteOrder.BIG_ENDIAN);
                bb.putLong(totalSequence);
                bb.order(ByteOrder.LITTLE_ENDIAN);
            }

            for (int color = 0; color < numColors; color++) {
                String sampleName = record.getParentGraph().getColor(color).getSampleName();

                bb.putInt(sampleName.length());
                bb.put(sampleName.getBytes());
            }

            for (int color = 0; color < numColors; color++) {
                byte[] errorRate = new byte[16];
                bb.put(errorRate);
            }

            for (int color = 0; color < numColors; color++) {
                bb.put((byte) (record.getParentGraph().getColor(color).isTipClippingApplied() ? 1 : 0));
                bb.put((byte) (record.getParentGraph().getColor(color).isLowCovgSupernodesRemoved() ? 1 : 0));
                bb.put((byte) (record.getParentGraph().getColor(color).isLowCovgKmersRemoved() ? 1 : 0));
                bb.put((byte) (record.getParentGraph().getColor(color).isCleanedAgainstGraph() ? 1 : 0));
                bb.putInt(record.getParentGraph().getColor(color).getLowCovSupernodesThreshold());
                bb.putInt(record.getParentGraph().getColor(color).getLowCovKmerThreshold());

                String cleanedAgainst = record.getParentGraph().getColor(color).getCleanedAgainstGraphName();

                bb.putInt(cleanedAgainst.length());
                bb.put(cleanedAgainst.getBytes());
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
            initialize(record);
        }

        ByteBuffer bb = ByteBuffer.allocateDirect(8 * this.numColors*4 + this.numColors);
        bb.clear();

        bb.order(ByteOrder.BIG_ENDIAN);

        long[] binaryKmer = record.getKmer();
        for (int i = 0; i < this.kmerBits; i++) {
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
