package uk.ac.ox.well.indiana.utils.io.cortex.graph;

import com.carrotsearch.sizeof.RamUsageEstimator;
import it.unimi.dsi.io.ByteBufferInputStream;
import org.apache.commons.collections.map.LRUMap;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.utils.BinaryFile;
import uk.ac.ox.well.indiana.utils.io.utils.BinaryUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class CortexGraph implements Iterable<CortexRecord>, Iterator<CortexRecord> {
    private File cortexFile;
    private BinaryFile in;

    //private CortexIndex index;
    private CortexHeader header;

    private long recordSize;
    private long numRecords;
    private long dataOffset;
    private long recordsSeen = 0;

    private ByteBufferInputStream mappedRecordBuffer = null;
    private CortexRecord nextRecord = null;

    private LRUMap cache = null;

    public CortexGraph(String cortexFilePath) {
        this.cortexFile = new File(cortexFilePath);
        loadCortexGraph(this.cortexFile, false);
    }

    public CortexGraph(File cortexFile) {
        this.cortexFile = cortexFile;
        loadCortexGraph(this.cortexFile, false);
    }

    public CortexGraph(String cortexFilePath, boolean ignoreIndex) {
        this.cortexFile = new File(cortexFilePath);
        loadCortexGraph(this.cortexFile, ignoreIndex);
    }

    public CortexGraph(File cortexFile, boolean ignoreIndex) {
        this.cortexFile = cortexFile;
        loadCortexGraph(this.cortexFile, ignoreIndex);
    }

    private byte[] fixStringsWithEarlyTerminators(byte[] string) {
        // Sometimes the names have an early terminator character (a bug in the old CORTEX output format).
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

    private void loadCortexGraph(File cortexFile, boolean ignoreIndex) {
        try {
            in = new BinaryFile(cortexFile, "r");

            byte[] headerStart = new byte[6];
            in.read(headerStart);
            String headerStartStr = new String(headerStart);

            if (!headerStartStr.equalsIgnoreCase("CORTEX")) {
                throw new RuntimeException("The file '" + cortexFile.getAbsolutePath() + "' does not appear to be a Cortex graph");
            }

            header = new CortexHeader();

            header.setVersion(in.readUnsignedInt());

            if (header.getVersion() != 6) {
                throw new RuntimeException("The file '" + cortexFile.getAbsolutePath() + "' is not a version 6 Cortex graph");
            }

            header.setKmerSize(in.readUnsignedInt());
            header.setKmerBits(in.readUnsignedInt());
            header.setNumColors(in.readUnsignedInt());

            for (int color = 0; color < header.getNumColors(); color++) {
                header.addColor(new CortexColor());
            }

            for (int color = 0; color < header.getNumColors(); color++) {
                header.getColor(color).setMeanReadLength(in.readUnsignedInt());
            }

            for (int color = 0; color < header.getNumColors(); color++) {
                header.getColor(color).setTotalSequence(in.readUnsignedLong());
            }

            for (int color = 0; color < header.getNumColors(); color++) {
                int sampleNameLength = in.readUnsignedInt();
                byte[] sampleName = new byte[sampleNameLength];
                in.read(sampleName);

                sampleName = fixStringsWithEarlyTerminators(sampleName);
                String sampleNameStr = new String(sampleName);

                header.getColor(color).setSampleName(sampleNameStr);
            }

            // Todo: fix this at some point - we're not actually getting the error rate properly
            for (int color = 0; color < header.getNumColors(); color++) {
                byte[] errorRate = new byte[16];
                in.read(errorRate);
            }

            for (int color = 0; color < header.getNumColors(); color++) {
                header.getColor(color).setTipClippingApplied(in.readBoolean());
                header.getColor(color).setLowCovgSupernodesRemoved(in.readBoolean());
                header.getColor(color).setLowCovgKmersRemoved(in.readBoolean());
                header.getColor(color).setCleanedAgainstGraph(in.readBoolean());
                header.getColor(color).setLowCovSupernodesThreshold(in.readUnsignedInt());
                header.getColor(color).setLowCovKmerThreshold(in.readUnsignedInt());

                int graphNameLength = in.readUnsignedInt();
                byte[] graphName = new byte[graphNameLength];
                in.read(graphName);

                graphName = fixStringsWithEarlyTerminators(graphName);

                header.getColor(color).setCleanedAgainstGraphName(new String(graphName));
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

            recordSize = (8*header.getKmerBits() + 5*header.getNumColors());
            numRecords = (dataSize / recordSize);

            mappedRecordBuffer = ByteBufferInputStream.map(in.getChannel(), FileChannel.MapMode.READ_ONLY);

            moveToBeginningOfRecordsSection();

            //index = new CortexIndex(ignoreIndex ? null : cortexFile.getAbsolutePath() + ".idx", this);

            cache = new LRUMap(100000);
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

    public long numRecordsSeen() { return recordsSeen; }

    private CortexRecord getNextRecord() {
        if (recordsSeen < getNumRecords()) {
            try {
                long[] binaryKmer = new long[header.getKmerBits()];
                for (int bits = 0; bits < header.getKmerBits(); bits++) {
                    byte[] b = new byte[8];

                    mappedRecordBuffer.read(b);

                    ByteBuffer bb = ByteBuffer.wrap(b);
                    binaryKmer[bits] = bb.getLong();
                }

                int[] coverages = new int[header.getNumColors()];
                for (int color = 0; color < header.getNumColors(); color++) {
                    byte[] coverage = new byte[4];
                    mappedRecordBuffer.read(coverage);

                    coverages[color] = BinaryUtils.toUnsignedInt(coverage);
                }

                byte[] edges = new byte[header.getNumColors()];
                mappedRecordBuffer.read(edges);

                recordsSeen++;

                return new CortexRecord(binaryKmer, coverages, edges, header.getKmerSize(), header.getKmerBits());
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

    public CortexRecord getRecord(long i) {
        long offset = dataOffset + i*recordSize;
        mappedRecordBuffer.position(offset);
        recordsSeen = 0;
        //return getNextRecord();

        CortexRecord cr = getNextRecord();

        if (cr != null) {
            cache.put(new ByteKmer(cr.getKmerAsBytes()), cr);
        }

        return cr;
    }

    public CortexRecord findRecord(ByteKmer kmer) {
        if (cache.containsKey(kmer)) {
            return (CortexRecord) cache.get(kmer);
        }

        long startIndex = 0;
        long stopIndex = getNumRecords();
        long midIndex = startIndex + (stopIndex - startIndex) / 2;

        while (startIndex != midIndex && midIndex != stopIndex) {
            long oldRecordsSeen = recordsSeen;

            CortexRecord startRecord = getRecord(startIndex);
            CortexRecord midRecord = getRecord(midIndex);
            CortexRecord stopRecord = getRecord(stopIndex);

            recordsSeen = oldRecordsSeen;

            ByteKmer startKmer = startRecord.getKmerAsByteKmer();
            ByteKmer midKmer = midRecord.getKmerAsByteKmer();
            ByteKmer stopKmer = stopRecord.getKmerAsByteKmer();

            if (startKmer.compareTo(stopKmer) > 0) {
                throw new IndianaException("Records are not sorted ('" + startKmer + "' is found before '" + stopKmer + "' but is lexicographically greater)");
            }

            if (startKmer.compareTo(midKmer) > 0) {
                throw new IndianaException("Records are not sorted ('" + startKmer + "' is found before '" + midKmer + "' but is lexicographically greater)");
            }

            if (kmer.compareTo(stopKmer) > 0 || kmer.compareTo(startKmer) < 0) { return null; }
            else if (startKmer.equals(kmer)) { return startRecord; }
            else if (midKmer.equals(kmer)) { return midRecord; }
            else if (stopKmer.equals(kmer)) { return stopRecord; }
            else if (kmer.compareTo(startKmer) > 0 && kmer.compareTo(midKmer) < 0) {
                stopIndex = midIndex;
                midIndex = startIndex + (stopIndex - startIndex) / 2;
            } else if (kmer.compareTo(midKmer) > 0 && kmer.compareTo(stopKmer) < 0) {
                startIndex = midIndex;
                midIndex = startIndex + ((stopIndex - startIndex) / 2);
            }
        }

        return null;
    }

    public CortexRecord findRecord(byte[] bk) { return findRecord(new ByteKmer(SequenceUtils.reverseComplement(bk))); }
    public CortexRecord findRecord(CortexKmer ck) { return findRecord(ck.getKmerAsBytes()); }
    public CortexRecord findRecord(String sk) { return findRecord(sk.getBytes()); }

    public File getCortexFile() { return cortexFile; }
    public CortexHeader getHeader() { return header; }
    public int getVersion() { return header.getVersion(); }
    public int getKmerSize() { return header.getKmerSize(); }
    public int getKmerBits() { return header.getKmerBits(); }
    public String getSampleName(int color) { return getColor(color).getSampleName(); }
    public int getNumColors() { return header.getNumColors(); }
    public long getNumRecords() { return numRecords; }
    public List<CortexColor> getColors() { return header.getColors(); }
    public boolean hasColor(int color) { return header.hasColor(color); }
    public CortexColor getColor(int color) { return header.getColor(color); }

    public int getColorForSampleName(String sampleName) {
        int sampleColor = -1;
        int sampleCopies = 0;

        for (int color = 0; color < header.getNumColors(); color++) {
            if (header.getColor(color).getSampleName().equalsIgnoreCase(sampleName)) {
                sampleColor = color;
                sampleCopies++;
            }
        }

        if (sampleColor == -1) {
            try {
                sampleColor = Integer.valueOf(sampleName);
                sampleCopies = 1;
            } catch (NumberFormatException e) {}
        }

        return (sampleCopies == 1) ? sampleColor : -1;
    }

    public List<Integer> getColorsForSampleNames(List<String> sampleNames) {
        List<Integer> colors = new ArrayList<>();

        if (sampleNames != null && !sampleNames.isEmpty()) {
            for (String sampleName : sampleNames) {
                colors.add(getColorForSampleName(sampleName));
            }
        }

        return colors;
    }

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
                    +  "  tip clipping: " + (cortexColor.isTipClippingApplied() ? "yes" : "no") + "\n"
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

}
