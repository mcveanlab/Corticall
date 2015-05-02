package uk.ac.ox.well.indiana.utils.io.kmerindex;

import ch.qos.logback.classic.Logger;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.indiana.Indiana;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.utils.BinaryFile;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

public class KmerIndex {
    public class KmerIndexRecord {
        private String kmer;
        private String chr;
        private int pos;

        public KmerIndexRecord(byte[] kmer, byte[] chr, int pos) {
            int chrNameLength = 0;
            for (int i = 0; i < chr.length; i++) {
                if (chr[i] == 0) {
                    chrNameLength = i;
                    break;
                }
            }

            byte[] chrName = new byte[chrNameLength];
            System.arraycopy(chr, 0, chrName, 0, chrNameLength);

            this.kmer = new String(kmer);
            this.chr = new String(chrName);
            this.pos = pos;
        }

        public String getKmer() { return kmer; }
        public String getChr() { return chr; }
        public int getPos() { return pos; }

        public String toString() { return String.format("%s %s:%d", kmer, chr, pos); }
    }

    private File file;
    private BinaryFile in;
    private final int bytesInChrName = 16;
    private final int bytesInPosField = 4;

    private int numRecords;
    private int kmerSize;
    private long dataOffset;
    private int recordSize;

    public KmerIndex(File file) {
        this.file = file;

        try {
            this.in = new BinaryFile(this.file, "r");

            byte[] headerStart = new byte[6];
            in.read(headerStart);
            String headerStartStr = new String(headerStart);

            if (!headerStartStr.equals("KINDEX")) {
                throw new IndianaException("The file '" + file.getAbsolutePath() + "' does not appear to be a valid kmer index");
            }

            numRecords = in.readUnsignedInt();
            kmerSize = in.readUnsignedInt();

            dataOffset = in.getFilePointer();
            recordSize = kmerSize + bytesInChrName + bytesInPosField;
        } catch (IOException e) {
            throw new IndianaException("IOException: ", e);
        }
    }

    public int getNumRecords() { return numRecords; }

    public KmerIndexRecord getRecord(int recordIndex) {
        long recordOffset = dataOffset + (recordIndex*recordSize);

        try {
            in.seek(recordOffset);

            byte[] kmer = new byte[kmerSize];
            byte[] chr = new byte[bytesInChrName];

            in.read(kmer);
            in.read(chr);
            int pos = in.readUnsignedInt();

            return new KmerIndexRecord(kmer, chr, pos);
        } catch (IOException e) {
            throw new IndianaException("Could not seek within index: ", e);
        }
    }

    public KmerIndexRecord findRecord(String kmer) {
        int startIndex = 0;
        int stopIndex = getNumRecords() - 1;
        int midIndex = startIndex + (stopIndex - startIndex) / 2;

        while (startIndex != midIndex && midIndex != stopIndex) {
            KmerIndexRecord startRecord = getRecord(startIndex);
            KmerIndexRecord midRecord   = getRecord(midIndex);
            KmerIndexRecord stopRecord  = getRecord(stopIndex);

            String startKmer = startRecord.getKmer();
            String midKmer   = midRecord.getKmer();
            String stopKmer  = stopRecord.getKmer();

            if (kmer.compareTo(stopKmer) > 0 || kmer.compareTo(startKmer) < 0) { return null; }
            else if (startKmer.equals(kmer)) { return startRecord; }
            else if (midKmer.equals(kmer))   { return midRecord; }
            else if (stopKmer.equals(kmer))  { return stopRecord; }
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

    public KmerIndexRecord findRecord(CortexKmer ck) {
        return findRecord(ck.getKmerAsString());
    }

    public Interval findLocus(String sk) {
        return findLocus(new CortexKmer(sk));
    }

    public Interval findLocus(CortexKmer ck) {
        KmerIndexRecord kir = findRecord(ck);

        return (kir == null) ? null : new Interval(kir.getChr(), kir.getPos(), kir.getPos() + kir.getKmer().length());
    }

    public static int writeIndex(File refFile, int kmerSize) {
        class Locus {
            private String chr;
            private int pos;

            public Locus(String chr, int pos) {
                this.chr = chr;
                this.pos = pos;
            }

            public String toString() {
                return String.format("%s:%d", chr, pos);
            }
        }

        Map<CortexKmer, Set<Locus>> kmerLoci = new TreeMap<CortexKmer, Set<Locus>>();

        File indexOut = new File(refFile.getAbsolutePath() + ".k" + kmerSize + ".idx");

        FastaSequenceFile ref = new FastaSequenceFile(refFile, true);
        ReferenceSequence rseq;
        while ((rseq = ref.nextSequence()) != null) {
            String[] name = rseq.getName().split("\\s+");
            String seq = new String(rseq.getBases());

            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                CortexKmer ck = new CortexKmer(seq.substring(i, i + kmerSize));
                Locus locus = new Locus(name[0], i);

                ContainerUtils.add(kmerLoci, ck, locus);
            }
        }

        int numKmers = 0;
        for (CortexKmer ck : kmerLoci.keySet()) {
            if (kmerLoci.get(ck).size() == 1) {
                numKmers++;
            }
        }

        try {
            FileOutputStream fos = new FileOutputStream(indexOut);
            FileChannel channel = fos.getChannel();

            ByteBuffer bb = ByteBuffer.allocateDirect(6 + 4 + 4);
            bb.order(ByteOrder.LITTLE_ENDIAN);
            bb.clear();

            bb.put("KINDEX".getBytes()); // magic number
            bb.putInt(numKmers);         // number of records in file
            bb.putInt(kmerSize);         // kmer size (number of bytes we expect to see for the kmer)

            bb.flip();

            channel.write(bb);

            bb = ByteBuffer.allocateDirect(numKmers * (kmerSize + 16 + 4));
            bb.order(ByteOrder.LITTLE_ENDIAN);
            bb.clear();

            int numKmersWritten = 0;
            for (CortexKmer ck : kmerLoci.keySet()) {
                if (kmerLoci.get(ck).size() == 1) {
                    Locus l = kmerLoci.get(ck).iterator().next();

                    byte[] buffer = new byte[16];
                    byte[] chr = l.chr.getBytes();

                    System.arraycopy(chr, 0, buffer, 0, chr.length);

                    bb.put(ck.getKmerAsBytes());
                    bb.put(buffer);
                    bb.putInt(l.pos);

                    numKmersWritten++;
                }
            }

            bb.flip();
            channel.write(bb);

            channel.close();
            fos.close();

            return numKmersWritten;
        } catch (FileNotFoundException e) {
            throw new IndianaException("File not found: " + e);
        } catch (IOException e) {
            throw new IndianaException("IOException: " + e);
        }
    }
}
