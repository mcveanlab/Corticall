package uk.ac.ox.well.indiana.commands.playground;

import htsjdk.samtools.*;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.Serializer;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.performance.PerformanceUtils;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.ConcurrentMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class index extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM")
    public File SAMFILE;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Output
    public PrintStream out;

    private long[] toLongArray(List<Chunk> chunks) {
        long[] la = new long[2*chunks.size()];

        for (int i = 0, j = 0; i < chunks.size(); i++, j += 2) {
            la[j]   = chunks.get(i).getChunkStart();
            la[j+1] = chunks.get(i).getChunkEnd();
        }

        return la;
    }

    private List<Chunk> toChunks(long[] la) {
        List<Chunk> chunks = new ArrayList<>();

        if (la != null && la.length >= 2) {
            for (int i = 0; i < la.length; i += 2) {
                Chunk chunk = new Chunk(la[i], la[i + 1]);

                chunks.add(chunk);
            }
        }

        return chunks;
    }

    @Override
    public void execute() {
        File kindex = new File(SAMFILE.getAbsolutePath().replaceAll(".bam", ".k" + KMER_SIZE + ".kindex"));
        kindex.delete();

        log.info("kindex: {}", kindex.getAbsolutePath());

        //DB db = DBMaker.fileDB(kindex).make();
        DB db = DBMaker.memoryDB().make();

        ConcurrentMap<Long, long[]> map = db.hashMap("map")
                .keySerializer(Serializer.LONG)
                .valueSerializer(Serializer.LONG_ARRAY)
                .create();

        SamReader sreader = SamReaderFactory.make()
                .setOption(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, true)
                .setOption(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES, true)
                .open(SAMFILE);

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing reads...")
                .updateRecord(100000)
                .message("reads processed")
                .make(log);

        Map<long[], long[]> m = new HashMap<>();

        for (SAMRecord sr : sreader) {
            BAMFileSpan sfs = (BAMFileSpan) sr.getFileSource().getFilePointer();

            //log.info("sfs: {}", sfs);

            byte[] rs = sr.getReadBases();

            for (int i = 0; i <= rs.length - KMER_SIZE; i++) {
                byte[] kmer = new byte[KMER_SIZE];

                System.arraycopy(rs, i, kmer, 0, KMER_SIZE);

                kmer = SequenceUtils.alphanumericallyLowestOrientation(kmer);

                List<Chunk> nchunks = sfs.getChunks();
                List<Chunk> ochunks = toChunks(map.get(kmer));

                List<Chunk> achunks = new ArrayList<>();
                achunks.addAll(nchunks);
                achunks.addAll(ochunks);

                achunks = Chunk.optimizeChunkList(achunks, 0);

                //map.put(kmer, toLongArray(achunks));
                m.put(CortexUtils.encodeBinaryKmer(kmer), toLongArray(achunks));
            }

            pm.update("reads processed, map size " + map.size() + ", " + PerformanceUtils.getCompactMemoryUsageStats());
        }

        // How many keys with multiple
        //for (byte[])

        /*
        log.info("Final map size: {}", map.size());

        log.info("Found {} reads", readsContainingKmer.size());

        chunks = Chunk.optimizeChunkList(chunks, 0);

        try {
            sreader.close();
        } catch (IOException e) {
            throw new IndianaException("Couldn't close file", e);
        }

        sreader = SamReaderFactory.makeDefault()
                .open(SAMFILE);

        int numReads = 0;
        for (Chunk chunk : chunks) {
            log.info("chunk: {}", chunk);
            SAMFileSpan sfs = new BAMFileSpan(chunk);

            SAMRecordIterator sri = sreader.indexing().iterator(sfs);

            while (sri.hasNext()) {
                SAMRecord sr = sri.next();

                numReads++;
            }

            log.info("Recovered {} reads from chunks", numReads);
        }

        log.info("");
        */
    }
}
