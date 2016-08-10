package uk.ac.ox.well.indiana.commands.playground;

import htsjdk.samtools.*;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexBinaryKmer;
import uk.ac.ox.well.indiana.utils.performance.PerformanceUtils;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class BuildIndex extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM")
    public File SAM_FILE;

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
        SamReader sreader = SamReaderFactory.make()
            .setOption(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, true)
            .setOption(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES, true)
            .open(SAM_FILE);

        KmerIndex ki = new KmerIndex(SAM_FILE, KMER_SIZE, true);

        ProgressMeter pm = new ProgressMeterFactory()
            .header("Processing reads...")
            .updateRecord(100000)
            .message("reads processed")
            .make(log);

        Map<CortexBinaryKmer, long[]> m = new TreeMap<>();
        long numKmersSeen = 0;

        for (SAMRecord sr : sreader) {
            BAMFileSpan sfs = (BAMFileSpan) sr.getFileSource().getFilePointer();

            byte[] rs = sr.getReadBases();

            for (int i = 0; i <= rs.length - KMER_SIZE; i++) {
                byte[] kmer = new byte[KMER_SIZE];

                System.arraycopy(rs, i, kmer, 0, KMER_SIZE);

                boolean kmerHasNs = false;
                for (int j = 0; j < kmer.length; j++) {
                    if (kmer[j] == 'N' || kmer[j] == 'n' || kmer[j] == '.') {
                        kmerHasNs = true;
                    }
                }

                if (!kmerHasNs) {
                    numKmersSeen++;

                    CortexBinaryKmer bk = new CortexBinaryKmer(kmer);

                    List<Chunk> nchunks = sfs.getChunks();

                    if (!m.containsKey(bk)) {
                        m.put(bk, toLongArray(nchunks));
                    } else {
                        List<Chunk> ochunks = toChunks(m.get(bk));

                        List<Chunk> achunks = new ArrayList<>();
                        achunks.addAll(nchunks);
                        achunks.addAll(ochunks);

                        achunks = Chunk.optimizeChunkList(achunks, 0);

                        m.put(bk, toLongArray(achunks));
                    }
                }
            }

            pm.update("reads processed, kmers seen " + numKmersSeen + ", map size " + m.size() + ", " + PerformanceUtils.getCompactMemoryUsageStats());
        }

        log.info("Writing index...");

        ki.putAll(m);
    }
}
