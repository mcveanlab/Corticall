package uk.ac.ox.well.cortexjdk.commands.assemble;

import com.carrotsearch.sizeof.RamUsageEstimator;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.*;
import uk.ac.ox.well.cortexjdk.utils.math.MoreMathUtils;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.*;

public class Build extends Module {
    @Argument(fullName="reads", shortName="r", doc="Read sequences")
    public ArrayList<File> READ_FILES;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE;

    @Argument(fullName="sampleName", shortName="s", doc="Sample name")
    public String SAMPLE_NAME;

    @Argument(fullName="threads", shortName="t", doc="Number of threads")
    public Integer NUM_THREADS = 1;

    @Output
    public File out;

    @Override
    public void execute() {
        verifyFiles(READ_FILES);

        log.info("Loading and sorting graph chunks...");
        Pair<List<CortexGraph>, Pair<Long, Integer>> val = buildSortedTempGraphs(READ_FILES, KMER_SIZE, SAMPLE_NAME, NUM_THREADS, out);

        log.info("Merging sorted graph chunks...");
        mergeSortedTempGraphs(val.getFirst(), val.getSecond());
    }

    private void verifyFiles(List<File> readFiles) {
        for (File readFile : readFiles) {
            String[] pieces = readFile.getAbsolutePath().split(":");

            for (String piece : pieces) {
                if (! new File(piece).exists()) {
                    throw new CortexJDKException("File not found: '" + piece + "'");
                }
            }
        }
    }

    private long presumableFreeMemory() {
        long allocatedMemory = (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory());

        return Runtime.getRuntime().maxMemory() - allocatedMemory;
    }

    private Pair<List<CortexGraph>, Pair<Long, Integer>> buildSortedTempGraphs(List<File> readFiles, int kmerSize, String sampleName, int numThreads, File cgout) {
        final long assumedRecordSize = RamUsageEstimator.sizeOf(new CortexRecord(CortexRecord.encodeBinaryKmer(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(kmerSize)),
                                                       new int[1],
                                                       new byte[1],
                                                       kmerSize,
                                                       CortexRecord.getKmerBits(kmerSize)));

        long freeMem = presumableFreeMemory();
        int overheadFactor = 2;
        int maxRecordsPerThread = (int) ((float) freeMem / (float) (overheadFactor*assumedRecordSize*numThreads));

        log.info("  using {} threads, max {} records per thread", numThreads, maxRecordsPerThread);

        ThreadPoolExecutor execIO   = (ThreadPoolExecutor) Executors.newFixedThreadPool(1);
        ThreadPoolExecutor execSort = (ThreadPoolExecutor) Executors.newFixedThreadPool(numThreads);

        BlockingQueue<List<CortexRecord>> readQueue = new ArrayBlockingQueue<>(2);
        BlockingQueue<List<CortexRecord>> writeQueue = new ArrayBlockingQueue<>(2);

        Future<Pair<Long, Integer>> gStatsFuture = execIO.submit(new RecordsProducer(readQueue, readFiles, kmerSize, maxRecordsPerThread, log));

        List<Future<List<CortexGraph>>> futures = new ArrayList<>();
        futures.add(execIO.submit(new RecordsWriter(out, sampleName, kmerSize, writeQueue, numThreads, log)));

        for (int i = 0; i < numThreads; i++) {
            execSort.submit(new RecordsSorter(readQueue, writeQueue, log));
        }

        List<CortexGraph> pieces = new ArrayList<>();
        Pair<Long, Integer> gstats;
        try {
            for (Future<List<CortexGraph>> future : futures) {
                pieces.addAll(future.get());
            }
            gstats = gStatsFuture.get();
        } catch (InterruptedException | ExecutionException e) {
            throw new CortexJDKException("Error retrieving temp graph results", e);
        }

        execIO.shutdown();
        execSort.shutdown();

        return new Pair<>(pieces, gstats);
    }

    private CortexGraph mergeSortedTempGraphs(List<CortexGraph> pieces, Pair<Long, Integer> gstats) {
        log.info("  merging {} chunks", pieces.size());

        CortexGraphWriter cgw = new CortexGraphWriter(out);

        CortexHeader ch = pieces.get(0).getHeader();
        ch.getColor(0).setTotalSequence(gstats.getFirst());
        ch.getColor(0).setMeanReadLength(gstats.getSecond());
        ch.getColor(0).setErrorRate(0.01);
        cgw.setHeader(ch);

        int kmerBits = CortexRecord.getKmerBits(KMER_SIZE);

        CortexCollection cc = new CortexCollection(pieces);

        long[] bk = null;
        int[] cov = new int[1];
        byte[] edges = new byte[1];
        for (CortexRecord cr : cc) {
            if (bk == null || !Arrays.equals(bk, cr.getBinaryKmer())) {
                if (bk != null) {
                    cgw.addRecord(new CortexRecord(bk, cov, edges, KMER_SIZE, kmerBits));
                }

                bk = cr.getBinaryKmer();
                cov[0] = 0;
                edges[0] = 0;

                cov[0] = MoreMathUtils.sum(cr.getCoverages());
                for (byte ed : cr.getEdges()) {
                    edges[0] |= ed;
                }
            } else if (Arrays.equals(bk, cr.getBinaryKmer())) {
                cov[0] += MoreMathUtils.sum(cr.getCoverages());
                for (byte ed : cr.getEdges()) {
                    edges[0] |= ed;
                }
            }
        }

        cgw.addRecord(new CortexRecord(bk, cov, edges, KMER_SIZE, kmerBits));
        cgw.close();

        CortexGraph cg = new CortexGraph(out);

        log.info("  -  wrote: {} unique records to {}", cg.getNumRecords(), cg.getFile().getAbsolutePath());

        return cg;
    }

}