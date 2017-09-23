package uk.ac.ox.well.cortexjdk.utils.alignment.kmer;

import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexBinaryKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.utils.BinaryFile;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by kiran on 20/09/2017.
 */
public class KmerStream {
    private BinaryFile bf;
    private int kmerSize;
    private int numBitFields;
    private String source;

    private int recordSize;
    private long dataOffset;
    private long numRecords;

    private KmerStreamRec nextRecord;
    private long nextRecordIndex;

    public KmerStream(File f) {
        try {
            this.bf = new BinaryFile(f, "r");

            byte[] magicWordStart = new byte[6];
            bf.read(magicWordStart);

            int version = bf.readInt();
            if (version != 1) {
                throw new CortexJDKException("Reference index has version " + version + ", expected version 1");
            }

            kmerSize = bf.readInt();

            int slength = bf.readInt();
            byte[] bsource = new byte[slength];
            bf.read(bsource);
            source = new String(bsource);

            byte[] magicWordEnd = new byte[6];
            bf.read(magicWordEnd);

            if (!Arrays.equals(magicWordStart, magicWordEnd)) {
                throw new CortexJDKException("Kmer index magic words mismatch");
            }

            this.numBitFields = CortexRecord.getKmerBits(kmerSize);
            this.recordSize = 8*numBitFields + 4 + 4;
            this.dataOffset = bf.getFilePointer();
            this.numRecords = (bf.length() - dataOffset) / recordSize;

            position(0);
        } catch (FileNotFoundException e) {
            throw new CortexJDKException("File not found", e);
        } catch (IOException e) {
            throw new CortexJDKException("IOException", e);
        }
    }

    public int getKmerSize() { return kmerSize; }

    public String getSource() { return source; }

    public long getNumRecords() {
        return numRecords;
    }

    private List<KmerStreamRec> adjacentRecords(long index) {
        KmerStreamRec ksr = position(index);

        List<KmerStreamRec> recs = new ArrayList<>();
        recs.add(ksr);

        for (long i = index - 1; i >= 0; i--) {
            KmerStreamRec ksri = position(i);

            if (ksr.kmerCompare(ksri) == 0) {
                recs.add(0, ksri);
            } else {
                break;
            }
        }

        for (long i = index + 1; i < numRecords - 1; i++) {
            KmerStreamRec ksri = position(i);

            if (ksr.kmerCompare(ksri) == 0) {
                recs.add(ksri);
            } else {
                break;
            }
        }

        return recs;
    }

    public List<KmerStreamRec> find(String kmer) {
        CortexBinaryKmer cbk = new CortexBinaryKmer(kmer.getBytes());

        long startIndex = 0;
        long stopIndex = numRecords - 1;
        long midIndex = startIndex + (stopIndex - startIndex) / 2;

        while (startIndex != midIndex && midIndex != stopIndex) {
            KmerStreamRec startRecord = position(startIndex);
            KmerStreamRec midRecord = position(midIndex);
            KmerStreamRec stopRecord = position(stopIndex);

            CortexBinaryKmer startKmer = startRecord.cbk;
            CortexBinaryKmer midKmer = midRecord.cbk;
            CortexBinaryKmer stopKmer = stopRecord.cbk;

            if (startKmer.compareTo(stopKmer) > 0) {
                throw new CortexJDKException("Records are not sorted ('" + startKmer + "' is found before '" + stopKmer + "' but is lexicographically greater)");
            }

            if (startKmer.compareTo(midKmer) > 0) {
                throw new CortexJDKException("Records are not sorted ('" + startKmer + "' is found before '" + midKmer + "' but is lexicographically greater)");
            }

            if (cbk.compareTo(stopKmer) > 0 || cbk.compareTo(startKmer) < 0) { return null; }
            else if (startKmer.equals(cbk)) { return adjacentRecords(startIndex); }
            else if (midKmer.equals(cbk)) { return adjacentRecords(midIndex); }
            else if (stopKmer.equals(cbk)) { return adjacentRecords(stopIndex); }
            else if (cbk.compareTo(startKmer) > 0 && cbk.compareTo(midKmer) < 0) {
                stopIndex = midIndex;
                midIndex = startIndex + (stopIndex - startIndex) / 2;
            } else if (cbk.compareTo(midKmer) > 0 && cbk.compareTo(stopKmer) < 0) {
                startIndex = midIndex;
                midIndex = startIndex + ((stopIndex - startIndex) / 2);
            }
        }

        return null;
    }

    private KmerStreamRec read() {
        try {
            long[] l = new long[numBitFields];
            for (int i = 0; i < l.length; i++) {
                l[i] = bf.readUnsignedLong();
            }
            int c = bf.readInt();
            int p = bf.readInt();

            nextRecordIndex++;

            return new KmerStreamRec(l, c, p);
        } catch (IOException e) {
            throw new CortexJDKException("Could not read: " + nextRecordIndex + " " + numRecords, e);
        }
    }

    private KmerStreamRec position(long i) {
        try {
            bf.seek(dataOffset + (i * recordSize));

            nextRecord = read();
            nextRecordIndex = i;

            return nextRecord;
        } catch (IOException e) {
            throw new CortexJDKException("Could not seek", e);
        }
    }

    public boolean hasNext() {
        return nextRecordIndex < numRecords - 1;
    }

    public KmerStreamRec peek() {
        return nextRecord;
    }

    public void advance() {
        nextRecord = hasNext() ? read() : null;
    }
}
