package uk.ac.ox.well.cortexjdk.commands.index.alignedbam;

import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexBinaryKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;

import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;

public class KmerIndex {
    private Path path;
    private int kmerSize;
    private int kmerBits;
    private FileChannel fc;

    private long dataOffset = 0;
    private int recordSize;
    private long numRecords;

    public KmerIndex(File bamFile, int kmerSize, boolean write) {
        String indexPath = bamFile.getAbsolutePath().replaceAll(".bam$", ".k" + kmerSize + "index");

        try {
            path = FileSystems.getDefault().getPath(indexPath);

            if (write) {
                if (new File(indexPath).exists()) {
                    throw new CortexJDKException("KmerIndex file '" + path + "' already exists.");
                }

                fc = FileChannel.open(path, EnumSet.of(StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE));

                writeHeader(kmerSize);
            } else {
                if (! new File(indexPath).exists()) {
                    throw new CortexJDKException("KmerIndex file '" + path + "' does not exist.");
                }

                fc = FileChannel.open(FileSystems.getDefault().getPath(indexPath), EnumSet.of(StandardOpenOption.READ));

                readHeader();
            }
        } catch (IOException e) {
            throw new CortexJDKException("Could not build path to KmerIndex file", e);
        }
    }

    private void writeHeader(int kmerSize) {
        try {
            this.kmerSize = kmerSize;
            this.kmerBits = CortexRecord.getKmerBits(kmerSize);
            this.recordSize = (8*kmerBits) + 8 + 8;

            ByteBuffer buffer = ByteBuffer.allocate(9 + 4 + 4);
            buffer.put("KMERINDEX".getBytes());
            buffer.putInt(this.kmerSize);
            buffer.putInt(this.kmerBits);

            buffer.flip();

            fc.write(buffer);

            //gos.write(buffer.array());
        } catch (IOException e) {
            throw new CortexJDKException("crapspasm", e);
        }
    }

    private void readHeader() {
        ByteBuffer magicNumberBuffer = ByteBuffer.allocate(9);
        ByteBuffer kmerSizeBuffer = ByteBuffer.allocate(4);
        ByteBuffer kmerBitsBuffer = ByteBuffer.allocate(4);

        try {
            fc.read(magicNumberBuffer);
            fc.read(kmerSizeBuffer);
            fc.read(kmerBitsBuffer);

            magicNumberBuffer.flip();
            kmerSizeBuffer.flip();
            kmerBitsBuffer.flip();

            String magicNumber = new String(magicNumberBuffer.array());

            if (!magicNumber.equals("KMERINDEX")) {
                throw new CortexJDKException("Not a kmer index file: '" + path + "'");
            }

            this.kmerSize = kmerSizeBuffer.asIntBuffer().get();
            this.kmerBits = kmerBitsBuffer.asIntBuffer().get();

            dataOffset = fc.position();

            recordSize = (8*kmerBits) + 8 + 8;

            numRecords = (fc.size() - dataOffset) / recordSize;
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public long getNumRecords() {
        return numRecords;
    }

    public void putAll(Map<CortexBinaryKmer, long[]> m) {
        int numEntries;
        int bufferSize;
        ByteBuffer buffer;

        try {
            for (CortexBinaryKmer bk : m.keySet()) {
                long[] locs = m.get(bk);

                numEntries = locs.length / 2;
                bufferSize = numEntries*recordSize;

                buffer = ByteBuffer.allocate(bufferSize);

                for (int i = 0; i < locs.length - 1; i += 2) {
                    for (int l = 0; l < kmerBits; l++) {
                        buffer.putLong(bk.getBinaryKmer()[l]);
                    }
                    buffer.putLong(locs[i]);
                    buffer.putLong(locs[i+1]);
                }

                buffer.flip();

                fc.write(buffer);
            }
        } catch (IOException e) {
            throw new CortexJDKException("crapspasm", e);
        }
    }

    public void close() {
        try {
            fc.close();
        } catch (IOException e) {
            throw new CortexJDKException("crapspasm", e);
        }
    }

    private long getFcPositionForIndex(long i) {
        return dataOffset + (i*recordSize);
    }

    public Pair<CortexBinaryKmer, long[]> getRecord(long i) {
        try {
            ByteBuffer recordBuffer = ByteBuffer.allocate(recordSize);

            fc.position(getFcPositionForIndex(i));

            fc.read(recordBuffer);
            long[] bk = new long[kmerBits];
            for (int b = 0; b < kmerBits; b++) {
                bk[b] = recordBuffer.getLong(8*b);
            }

            long[] l = new long[2];
            l[0] = recordBuffer.getLong(8*kmerBits);
            l[1] = recordBuffer.getLong((8*kmerBits) + (8));

            return new Pair<>(new CortexBinaryKmer(bk), l);
        } catch (IOException e) {
            throw new CortexJDKException("crapspasm", e);
        }
    }

    private List<long[]> makeList(long index) {
        Pair<CortexBinaryKmer, long[]> frecord = getRecord(index);

        List<long[]> chunks = new ArrayList<>();
        chunks.add(frecord.getSecond());

        long pindex = index;
        while (index > 0) {
            pindex--;

            Pair<CortexBinaryKmer, long[]> precord = getRecord(pindex);

            if (precord != null && precord.getFirst().equals(frecord.getFirst())) {
                chunks.add(0, precord.getSecond());
            } else {
                break;
            }
        }

        long nindex = index;
        while (index < numRecords) {
            nindex++;

            Pair<CortexBinaryKmer, long[]> nrecord = getRecord(nindex);

            if (nrecord != null && nrecord.getFirst().equals(frecord.getFirst())) {
                chunks.add(nrecord.getSecond());
            } else {
                break;
            }
        }

        return chunks;
    }

    public List<long[]> find(CortexBinaryKmer kmer) {
        long startIndex = 0;
        long stopIndex = numRecords - 1;
        long midIndex = startIndex + (stopIndex - startIndex) / 2;

        while (startIndex != midIndex && midIndex != stopIndex) {
            Pair<CortexBinaryKmer, long[]> startRecord = getRecord(startIndex);
            Pair<CortexBinaryKmer, long[]> midRecord = getRecord(midIndex);
            Pair<CortexBinaryKmer, long[]> stopRecord = getRecord(stopIndex);

            CortexBinaryKmer startKmer = startRecord.getFirst();
            CortexBinaryKmer midKmer   = midRecord.getFirst();
            CortexBinaryKmer stopKmer  = stopRecord.getFirst();

            if (kmer.compareTo(stopKmer) > 0 || kmer.compareTo(startKmer) < 0) {
                return null;
            } else if (startKmer.equals(kmer)) {
                return makeList(startIndex);
            } else if (midKmer.equals(kmer)) {
                return makeList(midIndex);
            } else if (stopKmer.equals(kmer)) {
                return makeList(stopIndex);
            } else if (kmer.compareTo(startKmer) > 0 && kmer.compareTo(midKmer) < 0) {
                stopIndex = midIndex;
                midIndex = startIndex + (stopIndex - startIndex) / 2;
            } else if (kmer.compareTo(midKmer) > 0 && kmer.compareTo(stopKmer) < 0) {
                startIndex = midIndex;
                midIndex = startIndex + ((stopIndex - startIndex) / 2);
            }
        }

        return null;
    }

    public List<long[]> find(byte[] kmer) {
        return find(new CortexBinaryKmer(kmer));
    }

    public List<long[]> find(CanonicalKmer ck) {
        return find(ck.getKmerAsBytes());
    }

    public List<long[]> find(String sk) { return find(sk.getBytes()); }
}
