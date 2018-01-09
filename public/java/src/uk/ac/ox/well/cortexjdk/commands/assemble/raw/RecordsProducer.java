package uk.ac.ox.well.cortexjdk.commands.assemble.raw;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.fastq.FastqRecord;
import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.reads.Reads;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;

import org.slf4j.Logger;
import uk.ac.ox.well.cortexjdk.utils.statistics.misc.StatisticsOnStream;

public class RecordsProducer implements Callable<Pair<Long, Integer>> {
    private BlockingQueue<List<CortexRecord>> queue;
    private List<File> readFiles;
    private int kmerSize;
    private int maxReads;
    private Logger log;

    public RecordsProducer(BlockingQueue<List<CortexRecord>> queue, List<File> readFiles, int kmerSize, int maxReads, Logger log) {
        this.readFiles = readFiles;
        this.kmerSize = kmerSize;
        this.maxReads = maxReads;
        this.queue = queue;
        this.log = log;
    }

    @Override
    public Pair<Long, Integer> call() throws Exception {
        List<CortexRecord> frs = new ArrayList<>(maxReads);

        try {
            StatisticsOnStream sos = new StatisticsOnStream();
            long totalSequence = 0;
            long numReads = 0;

            for (File readFile : readFiles) {
                Reads reads = new Reads(readFile);

                for (Pair<FastqRecord, FastqRecord> p : reads) {
                    List<String> rs = new ArrayList<>();
                    if (p.getFirst() != null) { rs.add(p.getFirst().getReadString()); }
                    if (p.getSecond() != null) { rs.add(p.getSecond().getReadString()); }

                    for (String r : rs) {
                        sos.push(r.length());
                        totalSequence += r.length();

                        for (int i = 0; i <= r.length() - kmerSize; i++) {
                            String sk = r.substring(i, i + kmerSize);

                            if (!sk.contains("N")) {
                                List<Set<String>> allInEdges = createSingleSampleInEdgeList(r, i);
                                List<Set<String>> allOutEdges = createSingleSampleOutEdgeList(r, i, kmerSize);

                                CortexRecord cr = new CortexRecord(sk, Lists.newArrayList(1), allInEdges, allOutEdges);

                                frs.add(cr);
                            }
                        }

                        numReads++;
                    }

                    if (frs.size() >= maxReads) {
                        queue.put(frs);

                        log.info("  - queued: {} records, {} reads processed [{}]", frs.size(), numReads, Thread.currentThread().getName());

                        frs = new ArrayList<>();
                    }
                }
            }

            if (frs.size() > 0) {
                queue.put(frs);

                log.info("  - queued: {} records, {} reads processed [{}]", frs.size(), numReads, Thread.currentThread().getName());
            }

            queue.put(new ArrayList<>());

            return new Pair<>(totalSequence, (int) Math.floor(sos.getMean()));
        } catch (InterruptedException e) {
            throw new CortexJDKException("Interrupted while adding records to processing queue", e);
        } catch (FileNotFoundException e) {
            throw new CortexJDKException("File not found", e);
        }
    }

    @NotNull
    private List<Set<String>> createSingleSampleOutEdgeList(String r, int i, int kmerSize) {
        Set<String> outEdges = i == r.length() - kmerSize ? Sets.newHashSet() : Sets.newHashSet(r.substring(i + kmerSize, i + kmerSize + 1));
        List<Set<String>> allOutEdges = new ArrayList<>();
        allOutEdges.add(outEdges);
        return allOutEdges;
    }

    @NotNull
    private List<Set<String>> createSingleSampleInEdgeList(String r, int i) {
        Set<String> inEdges = i == 0 ? Sets.newHashSet() : Sets.newHashSet(r.substring(i - 1, i));
        List<Set<String>> allInEdges = new ArrayList<>();
        allInEdges.add(inEdges);
        return allInEdges;
    }
}
