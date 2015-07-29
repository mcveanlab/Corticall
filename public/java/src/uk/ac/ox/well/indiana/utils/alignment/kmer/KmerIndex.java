package uk.ac.ox.well.indiana.utils.alignment.kmer;

import com.javacodegeeks.stringsearch.BM;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Interval;
import it.unimi.dsi.io.ByteBufferInputStream;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.utils.BinaryFile;
import uk.ac.ox.well.indiana.utils.io.utils.BinaryUtils;
import uk.ac.ox.well.indiana.utils.performance.PerformanceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.*;

public class KmerIndex {
    private File kIndexFile;
    private ByteBufferInputStream mappedRecordBuffer = null;

    private int version;
    private int kmerSize;
    private int numKmers;

    private int recordSize;
    private long numRecords;
    private long dataOffset;

    private Map<Byte, String> dictMap;
    private Map<Byte, String> refMap;

    private Map<String, KmerIndexRecord> kirCache = new HashMap<String, KmerIndexRecord>();

    public class KmerIndexRecord {
        private String sk;
        private Set<Interval> intervals = new HashSet<Interval>();
        private boolean isUnique = true;

        public KmerIndexRecord(byte[] kmer) {
            sk = new String(kmer);
        }

        public KmerIndexRecord(byte[] kmer, byte chr, int pos) {
            sk = new String(kmer);
            addInterval(chr, pos);
        }

        public void addInterval(byte chr, int pos) {
            if (chr == -1) {
                isUnique = false;
            } else {
                intervals.add(new Interval(dictMap.get(chr), pos, pos));
            }
        }

        public String getKmer() { return sk; }
        public Set<Interval> getLocations() { return intervals; }
        public boolean isUnique() { return isUnique; }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            KmerIndexRecord that = (KmerIndexRecord) o;

            if (isUnique != that.isUnique) return false;
            if (sk != null ? !sk.equals(that.sk) : that.sk != null) return false;
            return !(intervals != null ? !intervals.equals(that.intervals) : that.intervals != null);

        }

        @Override
        public int hashCode() {
            int result = sk != null ? sk.hashCode() : 0;
            result = 31 * result + (intervals != null ? intervals.hashCode() : 0);
            result = 31 * result + (isUnique ? 1 : 0);
            return result;
        }

        @Override
        public String toString() {
            return "KmerIndexRecord{" +
                    "sk='" + sk + '\'' +
                    ", intervals=" + intervals +
                    ", isUnique=" + isUnique +
                    '}';
        }
    }

    public KmerIndex(File refFile, int kmerSize) {
        this.kIndexFile = new File(refFile.getAbsoluteFile() + ".k" + kmerSize + "index");

        if (!kIndexFile.exists()) {
            createIndex(refFile, kmerSize, kIndexFile);
        }

        loadDictAndRef(refFile);

        loadIndex(refFile, kIndexFile);
    }

    private void loadDictAndRef(File refFile) {
        dictMap = new HashMap<Byte, String>();
        refMap = new HashMap<Byte, String>();

        FastaSequenceFile ref = new FastaSequenceFile(refFile, true);
        ReferenceSequence rseq;

        while ((rseq = ref.nextSequence()) != null) {
            String name = rseq.getName().split("\\s+")[0];
            byte chr = (byte) rseq.getContigIndex();

            dictMap.put(chr, name);
            refMap.put(chr, new String(rseq.getBases()));
        }
    }

    private void createIndex(File refFile, int kmerSize, File kIndexFile) {
        Map<String, GenomePosition> kmers = new TreeMap<String, GenomePosition>();

        ReferenceSequenceFile ref = new FastaSequenceFile(refFile, true);
        ReferenceSequence rseq;

        while ((rseq = ref.nextSequence()) != null) {
            String seq = new String(rseq.getBases());
            byte refIndex = (byte) rseq.getContigIndex();

            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                String sk = seq.substring(i, i + kmerSize);

                if (!kmers.containsKey(sk)) {
                    GenomePosition gp = new GenomePosition(refIndex, i);
                    kmers.put(sk, gp);
                } else if (kmers.containsKey(sk) && kmers.get(sk) != null) {
                    kmers.put(sk, null);
                }
            }
        }

        try {
            FileOutputStream fos = new FileOutputStream(kIndexFile);
            FileChannel channel = fos.getChannel();

            ByteBuffer bb = ByteBuffer.allocateDirect(6 + 4 + 4 + 4);
            bb.order(ByteOrder.LITTLE_ENDIAN);
            bb.clear();

            bb.put("KINDEX".getBytes());
            bb.putInt(0);
            bb.putInt(kmerSize);
            bb.putInt(kmers.size());

            bb.flip();
            channel.write(bb);

            bb = ByteBuffer.allocateDirect(kmerSize + 1 + 4);
            bb.order(ByteOrder.LITTLE_ENDIAN);

            for (String sk : kmers.keySet()) {
                bb.clear();

                bb.put(sk.getBytes());

                GenomePosition gp = kmers.get(sk);
                if (gp == null) {
                    bb.put((byte) -1);
                    bb.putInt(0);
                } else {
                    bb.put(gp.chr);
                    bb.putInt(gp.pos);
                }

                bb.flip();
                channel.write(bb);
            }

            channel.close();
            fos.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void loadIndex(File refFile, File kIndexFile) {
        try {
            BinaryFile in = new BinaryFile(kIndexFile, "r");

            byte[] headerStart = new byte[6];
            in.read(headerStart);
            String headerStartStr = new String(headerStart);

            this.version = in.readInt();
            this.kmerSize = in.readUnsignedInt();
            this.numKmers = in.readUnsignedInt();

            long size = in.getChannel().size();
            dataOffset = in.getFilePointer();
            long dataSize = size - dataOffset;

            this.recordSize = (kmerSize + 1 + 4);
            this.numRecords = (dataSize / recordSize);

            mappedRecordBuffer = ByteBufferInputStream.map(in.getChannel(), FileChannel.MapMode.READ_ONLY);
            moveToBeginningOfRecordsSection();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void moveToBeginningOfRecordsSection() {
        mappedRecordBuffer.position(dataOffset);
    }

    private KmerIndexRecord getNextRecord() {
        try {
            byte[] kmer = new byte[kmerSize];
            mappedRecordBuffer.read(kmer);

            byte chr = (byte) mappedRecordBuffer.read();

            byte[] posb = new byte[4];
            mappedRecordBuffer.read(posb);

            int pos = BinaryUtils.toUnsignedInt(posb);

            return new KmerIndexRecord(kmer, chr, pos);
        } catch (IOException e) {
            throw new IndianaException("Couldn't read next kmer index entry: ", e);
        }
    }

    private KmerIndexRecord getRecord(long i) {
        long offset = dataOffset + i*recordSize;
        mappedRecordBuffer.position(offset);
        return getNextRecord();
    }

    private KmerIndexRecord findRecord(byte[] bk) {
        return findRecord(new String(bk));
    }

    private KmerIndexRecord findRecord(String kmer) {
        long startIndex = 0;
        long stopIndex = numRecords - 1;
        long midIndex = startIndex + (stopIndex - startIndex) / 2;

        while (startIndex != midIndex && midIndex != stopIndex) {
            KmerIndexRecord startRecord = getRecord(startIndex);
            KmerIndexRecord midRecord = getRecord(midIndex);
            KmerIndexRecord stopRecord = getRecord(stopIndex);

            String startKmer = startRecord.getKmer();
            String midKmer = midRecord.getKmer();
            String stopKmer = stopRecord.getKmer();

            if (startKmer.compareTo(stopKmer) > 0) {
                throw new IndianaException("Records are not sorted ('" + startKmer + "' is found before '" + stopKmer + "' but is lexicographically greater)");
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

    private KmerIndexRecord slowLookup(byte[] kmer) {
        KmerIndexRecord kir = new KmerIndexRecord(kmer);

        String sk = kir.getKmer();

        for (Byte chr : refMap.keySet()) {
            String seq = refMap.get(chr);

            List<Integer> positions = BM.findAll(sk, seq);
            for (Integer pos : positions) {
                kir.addInterval(chr, pos);
            }
        }

        return kir.getLocations().size() > 0 ? kir : null;
    }

    public KmerIndexRecord lookup(byte[] kmer) {
        return lookup(new String(kmer));
    }

    public KmerIndexRecord lookup(String sk) {
        if (kirCache.containsKey(sk)) {
            return kirCache.get(sk);
        } else {
            KmerIndexRecord kir = findRecord(sk.getBytes());

            if (kir != null && !kir.isUnique()) {
                kir = slowLookup(sk.getBytes());

                if (kir != null) {
                    kirCache.put(kir.getKmer(), kir);
                }
            }

            return kir;
        }
    }
}
