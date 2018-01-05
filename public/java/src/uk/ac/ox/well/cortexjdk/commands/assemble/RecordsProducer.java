package uk.ac.ox.well.cortexjdk.commands.assemble;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.fastq.FastqRecord;
import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.reads.Reads;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.BlockingQueue;
import org.slf4j.Logger;

public class RecordsProducer implements Runnable {
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
    public void run() {
        List<CortexRecord> frs = new ArrayList<>(maxReads);

        try {
            for (File readFile : readFiles) {
                Reads reads = new Reads(readFile);

                for (Pair<FastqRecord, FastqRecord> p : reads) {
                    List<String> rs = new ArrayList<>();
                    if (p.getFirst() != null) { rs.add(p.getFirst().getReadString()); }
                    if (p.getSecond() != null) { rs.add(p.getSecond().getReadString()); }

                    for (String r : rs) {
                        for (int i = 0; i <= r.length() - kmerSize; i++) {
                            String sk = r.substring(i, i + kmerSize);

                            List<Set<String>> allInEdges = createSingleSampleInEdgeList(r, i);
                            List<Set<String>> allOutEdges = createSingleSampleOutEdgeList(r, i, kmerSize);

                            CortexRecord cr = new CortexRecord(sk, Lists.newArrayList(1), allInEdges, allOutEdges);

                            frs.add(cr);
                        }
                    }

                    if (frs.size() >= maxReads) {
                        queue.put(frs);

                        log.info("  queued {} records [{}]", frs.size(), Thread.currentThread().getName());

                        frs = new ArrayList<>();
                    }
                }
            }

            if (frs.size() > 0) {
                queue.put(frs);

                log.info("  queued last {} records [{}]", frs.size(), Thread.currentThread().getName());
            }

            queue.put(new ArrayList<>());
        } catch (InterruptedException e) {
            throw new CortexJDKException("Interrupted while adding records to processing queue", e);
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
