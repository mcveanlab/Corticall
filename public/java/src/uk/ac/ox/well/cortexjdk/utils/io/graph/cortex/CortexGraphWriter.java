package uk.ac.ox.well.cortexjdk.utils.io.graph.cortex;

import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;

public class CortexGraphWriter {
    private File cortexFile;
    private FileOutputStream fos = null;
    private FileChannel channel;

    private CortexHeader header;

    public CortexGraphWriter(File cortexFile) {
        this.cortexFile = cortexFile;
    }

    public CortexGraphWriter(String cortexFilePath) {
        this.cortexFile = new File(cortexFilePath);
    }

    public void setHeader(CortexHeader header) { this.header = header; }
    public CortexHeader getHeader() { return this.header; }

    private void initialize() {
        try {
            fos = new FileOutputStream(cortexFile);
            channel = fos.getChannel();

            int stringLengths = 0;
            for (CortexColor c : header.getColors()) {
                stringLengths += c.getSampleName().length() + c.getCleanedAgainstGraphName().length();
            }

            ByteBuffer bb = ByteBuffer.allocateDirect(2 * (44 + header.getNumColors()*32 + stringLengths));
            bb.order(ByteOrder.LITTLE_ENDIAN);
            bb.clear();

            bb.put("CORTEX".getBytes());

            bb.putInt(header.getVersion());
            bb.putInt(header.getKmerSize());
            bb.putInt(header.getKmerBits());
            bb.putInt(header.getNumColors());

            for (int color = 0; color < header.getNumColors(); color++) {
                int meanReadLength = header.getColor(color).getMeanReadLength();
                bb.putInt(meanReadLength);
            }

            for (int color = 0; color < header.getNumColors(); color++) {
                long totalSequence = header.getColor(color).getTotalSequence();

                bb.order(ByteOrder.BIG_ENDIAN);
                bb.putLong(totalSequence);
                bb.order(ByteOrder.LITTLE_ENDIAN);
            }

            for (int color = 0; color < header.getNumColors(); color++) {
                String sampleName = header.getColor(color).getSampleName();

                bb.putInt(sampleName.length());
                bb.put(sampleName.getBytes());
            }

            for (int color = 0; color < header.getNumColors(); color++) {
                byte[] errorRate = new byte[16];
                bb.put(errorRate);
            }

            for (int color = 0; color < header.getNumColors(); color++) {
                bb.put((byte) (header.getColor(color).isTipClippingApplied() ? 1 : 0));
                bb.put((byte) (header.getColor(color).isLowCovgSupernodesRemoved() ? 1 : 0));
                bb.put((byte) (header.getColor(color).isLowCovgKmersRemoved() ? 1 : 0));
                bb.put((byte) (header.getColor(color).isCleanedAgainstGraph() ? 1 : 0));
                bb.putInt(header.getColor(color).getLowCovSupernodesThreshold());
                bb.putInt(header.getColor(color).getLowCovKmerThreshold());

                String cleanedAgainst = header.getColor(color).getCleanedAgainstGraphName();

                bb.putInt(cleanedAgainst.length());
                bb.put(cleanedAgainst.getBytes());
            }

            bb.put("CORTEX".getBytes());

            bb.flip();

            channel.write(bb);
        } catch (FileNotFoundException e) {
            throw new CortexJDKException("Unable to open file '" + cortexFile.getAbsolutePath() + "'", e);
        } catch (IOException e) {
            throw new CortexJDKException("Unable to write header to file '" + cortexFile.getAbsolutePath() + "'", e);
        }
    }

    public void addRecord(CortexRecord record) {
        if (fos == null) { initialize(); }

        ByteBuffer bb = ByteBuffer.allocateDirect(8*header.getKmerBits() + header.getNumColors()*5);
        bb.clear();

        bb.order(ByteOrder.BIG_ENDIAN);

        long[] binaryKmer = record.getBinaryKmer();
        for (int i = 0; i < header.getKmerBits(); i++) {
            bb.putLong(binaryKmer[i]);
        }

        bb.order(ByteOrder.LITTLE_ENDIAN);

        int[] coverage = record.getCoverages();
        for (int c = 0; c < header.getNumColors(); c++) {
            bb.putInt(coverage[c]);
        }

        byte[] edges = record.getEdges();
        for (int c = 0; c < header.getNumColors(); c++) {
            bb.put(edges[c]);
        }

        bb.flip();

        try {
            channel.write(bb);
        } catch (IOException e) {
            throw new CortexJDKException("Unable to write record to file '" + cortexFile.getAbsolutePath() + "'", e);
        }
    }

    public void close() {
        if (fos == null) { initialize(); }

        try {
            channel.close();
            fos.close();
        } catch (IOException e) {
            throw new CortexJDKException("Unable to close '" + cortexFile.getAbsolutePath() + "'", e);
        }
    }
}
