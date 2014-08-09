package uk.ac.ox.well.indiana.utils.io.cortex.graph;

import com.carrotsearch.sizeof.RamUsageEstimator;
import it.unimi.dsi.io.ByteBufferInputStream;
import uk.ac.ox.well.indiana.utils.io.utils.BinaryFile;
import uk.ac.ox.well.indiana.utils.io.utils.BinaryUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;

public class CortexGraph implements Iterable<CortexRecord>, Iterator<CortexRecord> {
    private File cortexFile;
    private BinaryFile in;

    private int version;
    private int kmerSize;
    private int kmerBits;
    private int numColors;
    private long numRecords;

    private long dataOffset;
    private long recordSize;
    //private long sizeOfBuffer;
    private long recordsSeen = 0;

    private ArrayList<CortexColor> colors = new ArrayList<CortexColor>();
    //private HashMap<String, Integer> nameToColor = new HashMap<String, Integer>();

    private ByteBufferInputStream mappedRecordBuffer = null;
    private CortexRecord nextRecord = null;

    public CortexGraph(String cortexFilePath) {
        this.cortexFile = new File(cortexFilePath);
        loadCortexGraph(this.cortexFile);
    }

    public CortexGraph(File cortexFile) {
        this.cortexFile = cortexFile;
        loadCortexGraph(this.cortexFile);
    }

    private byte[] fixStringsWithEarlyTerminators(byte[] string) {
        // Sometimes the names have an early terminator character (a bug in the CORTEX output format).
        int earlyTerminatorPosition = string.length;
        for (int i = string.length - 1; i >= 0; i--) {
            if (string[i] == 0x0) {
                earlyTerminatorPosition = i;
            }
        }

        if (earlyTerminatorPosition < string.length) {
            return Arrays.copyOfRange(string, 0, earlyTerminatorPosition);
        }

        return string;
    }

    private void loadCortexGraph(File cortexFile) {
        try {
            in = new BinaryFile(cortexFile, "r");

            byte[] headerStart = new byte[6];
            in.read(headerStart);
            String headerStartStr = new String(headerStart);

            if (!headerStartStr.equalsIgnoreCase("CORTEX")) {
                throw new RuntimeException("The file '" + cortexFile.getAbsolutePath() + "' does not appear to be a Cortex graph");
            }

            version = in.readUnsignedInt();

            if (version != 6) {
                throw new RuntimeException("The file '" + cortexFile.getAbsolutePath() + "' is not a version 6 Cortex graph");
            }

            kmerSize = in.readUnsignedInt();
            kmerBits = in.readUnsignedInt();
            numColors = in.readUnsignedInt();

            for (int color = 0; color < numColors; color++) {
                colors.add(new CortexColor());
            }

            for (int color = 0; color < numColors; color++) {
                colors.get(color).setMeanReadLength(in.readUnsignedInt());
            }

            for (int color = 0; color < numColors; color++) {
                colors.get(color).setTotalSequence(in.readUnsignedLong());
            }

            for (int color = 0; color < numColors; color++) {
                int sampleNameLength = in.readUnsignedInt();
                byte[] sampleName = new byte[sampleNameLength];
                in.read(sampleName);

                sampleName = fixStringsWithEarlyTerminators(sampleName);
                String sampleNameStr = new String(sampleName);

                colors.get(color).setSampleName(sampleNameStr);
                //nameToColor.put(sampleNameStr, color);
            }

            // Todo: fix this at some point - we're not actually getting the error rate properly
            for (int color = 0; color < numColors; color++) {
                byte[] errorRate = new byte[16];
                in.read(errorRate);
            }

            for (int color = 0; color < numColors; color++) {
                colors.get(color).setTopClippingApplied(in.readBoolean());
                colors.get(color).setLowCovgSupernodesRemoved(in.readBoolean());
                colors.get(color).setLowCovgKmersRemoved(in.readBoolean());
                colors.get(color).setCleanedAgainstGraph(in.readBoolean());
                colors.get(color).setLowCovSupernodesThreshold(in.readUnsignedInt());
                colors.get(color).setLowCovKmerThreshold(in.readUnsignedInt());

                int graphNameLength = in.readUnsignedInt();
                byte[] graphName = new byte[graphNameLength];
                in.read(graphName);

                graphName = fixStringsWithEarlyTerminators(graphName);

                colors.get(color).setCleanedAgainstGraphName(new String(graphName));
            }

            byte[] headerEnd = new byte[6];
            in.read(headerEnd);
            String headerEndStr = new String(headerEnd);

            if (!headerEndStr.equalsIgnoreCase("CORTEX")) {
                throw new RuntimeException("We didn't see a proper header terminator at the expected place in Cortex graph '" + cortexFile.getAbsolutePath() + "'");
            }

            long size = in.getChannel().size();
            dataOffset = in.getFilePointer();
            long dataSize = size - dataOffset;

            recordSize = (8*kmerBits + 4*numColors + 1*numColors);
            numRecords = dataSize / recordSize;

            mappedRecordBuffer = ByteBufferInputStream.map(in.getChannel(), FileChannel.MapMode.READ_ONLY);
            //mappedRecordBuffer.position(dataOffset);

            //nextRecord = getNextRecord();
            moveToBeginningOfRecordsSection();
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Cortex graph file '" + cortexFile.getAbsolutePath() + "' not found: " + e);
        } catch (IOException e) {
            throw new RuntimeException("Error while parsing Cortex graph file '" + cortexFile.getAbsolutePath() + "': " + e);
        }
    }

    private void moveToBeginningOfRecordsSection() {
        mappedRecordBuffer.position(dataOffset);
        recordsSeen = 0;
        nextRecord = getNextRecord();
    }

    public File getCortexFile() { return cortexFile; }

    public int getVersion() { return version; }

    public int getKmerSize() { return kmerSize; }

    public int getKmerBits() { return kmerBits; }

    public int getNumColors() { return numColors; }

    public ArrayList<CortexColor> getColors() { return colors; }

    public String toString() {
        String info = "file: " + cortexFile.getAbsolutePath() + "\n"
                    + "----" + "\n"
                    + "binary version: " + this.getVersion() + "\n"
                    + "kmer size: " + this.getKmerSize() + "\n"
                    + "bitfields: " + this.getKmerBits() + "\n"
                    + "colors: " + this.getNumColors() + "\n";

        for (int color = 0; color < this.getNumColors(); color++) {
            CortexColor cortexColor = this.getColors().get(color);
            info += "-- Color " + color + " --\n"
                 +  "  sample name: '" + cortexColor.getSampleName() + "'\n"
                 +  "  mean read length: " + cortexColor.getMeanReadLength() + "\n"
                 +  "  total sequence loaded: " + "(not parsed)" + "\n"
                 +  "  sequence error rate: " + "(not parsed)" + "\n"
                 +  "  tip clipping: " + (cortexColor.isTopClippingApplied() ? "yes" : "no") + "\n"
                 +  "  remove_low_coverage_supernodes: " + (cortexColor.isLowCovgSupernodesRemoved() ? "yes" : "no") + "\n"
                 +  "  remove_low_coverage_kmers: " + (cortexColor.isLowCovgKmersRemoved() ? "yes" : "no") + "\n"
                 +  "  cleaned against graph: " + (cortexColor.isCleanedAgainstGraph() ? "yes" : "no") + "\n";
        }

        info += "----" + "\n";
        info += "kmers: " + getNumRecords() + "\n";
        info += "----" + "\n";
        info += "size of one Cortex record: " + RamUsageEstimator.humanSizeOf(nextRecord) + "\n";
        info += "size of full Cortex graph: " + RamUsageEstimator.humanReadableUnits(RamUsageEstimator.sizeOf(nextRecord)*getNumRecords()) + "\n";

        return info;
    }

    public long getNumRecords() {
        return numRecords;
    }

    private CortexRecord getNextRecord() {
        if (recordsSeen < numRecords) {
            try {
                long[] binaryKmer = new long[kmerBits];
                for (int bits = 0; bits < kmerBits; bits++) {
                    //binaryKmer[bits] = mappedRecordBuffer.getLong();

                    byte[] b = new byte[8];

                    mappedRecordBuffer.read(b);

                    ByteBuffer bb = ByteBuffer.wrap(b);
                    binaryKmer[bits] = bb.getLong();
                }

                int[] coverages = new int[numColors];
                for (int color = 0; color < numColors; color++) {
                    byte[] coverage = new byte[4];
                    //mappedRecordBuffer.get(coverage);
                    mappedRecordBuffer.read(coverage);

                    coverages[color] = BinaryUtils.toUnsignedInt(coverage);
                }

                byte[] edges = new byte[numColors];
                //for (int color = 0; color < numColors; color++) {
                    //edges[color] = mappedRecordBuffer.get();
                    //mappedRecordBuffer.read(edges[color]);
                //}
                mappedRecordBuffer.read(edges);

                recordsSeen++;

                return new CortexRecord(binaryKmer, coverages, edges, kmerSize, kmerBits);
            } catch (IOException e) {
                return null;
            }
        }

        return null;
    }

    public Iterator<CortexRecord> iterator() {
        moveToBeginningOfRecordsSection();

        return this;
    }

    public boolean hasNext() {
        return nextRecord != null;
    }

    public CortexRecord next() {
        CortexRecord currentRecord = nextRecord;

        nextRecord = getNextRecord();
        if (nextRecord == null) {
            close();
        }

        return currentRecord;
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }

    public void close() {
        try {
            this.in.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public boolean hasColor(int color) {
        return (color < colors.size());
    }

    public CortexColor getColor(int color) {
        return colors.get(color);
    }

    public int getColorForSampleName(String sampleName) {
        int sampleColor = -1;
        int sampleCopies = 0;

        for (int color = 0; color < colors.size(); color++) {
            if (colors.get(color).getSampleName().equalsIgnoreCase(sampleName)) {
                sampleColor = color;
                sampleCopies++;
            }
        }

        return (sampleCopies == 1) ? sampleColor : -1;
    }
}
