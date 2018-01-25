package uk.ac.ox.well.cortexjdk.commands.assemble.raw;

import org.jetbrains.annotations.NotNull;
import org.slf4j.Logger;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;

public class RecordsWriter implements Callable<List<CortexGraph>> {
    private File cgout;
    private String sampleName;
    private int kmerSize;
    private BlockingQueue<List<CortexRecord>> writeQueue;
    private int numThreads;
    private Logger log;

    private int numPoisonPills = 0;

    public RecordsWriter(File cgout, String sampleName, int kmerSize, BlockingQueue<List<CortexRecord>> writeQueue, int numThreads, Logger log) {
        this.cgout = cgout;
        this.sampleName = sampleName;
        this.kmerSize = kmerSize;
        this.writeQueue = writeQueue;
        this.numThreads = numThreads;
        this.log = log;
    }

    @NotNull
    private File makeTempGraphFile(File cgout) {
        File tempFile;
        try {
            tempFile = File.createTempFile("temp", ".ctx", cgout.getAbsoluteFile().getParentFile());
            tempFile.deleteOnExit();
        } catch (IOException e) {
            throw new CortexJDKException("Could not create temp file in directory '" + cgout.getAbsoluteFile().getParent() + "'");
        }
        return tempFile;
    }

    private CortexGraphWriter makeGraphWriter(String sampleName, int kmerSize, File cgout) {
        CortexColor cc = new CortexColor();
        cc.setSampleName(sampleName);

        CortexHeader ch = new CortexHeader();
        ch.setVersion(6);
        ch.setNumColors(1);
        ch.setKmerSize(kmerSize);
        ch.setKmerBits(CortexRecord.getKmerBits(kmerSize));
        ch.addColor(cc);

        CortexGraphWriter cgw = new CortexGraphWriter(cgout);
        cgw.setHeader(ch);

        return cgw;
    }

    @Override
    public List<CortexGraph> call() throws Exception {
        List<CortexGraph> graphs = new ArrayList<>();

        try {
            List<CortexRecord> lrs;

            while (true) {
                lrs = writeQueue.peek();

                if (numPoisonPills == numThreads) {
                    break;
                } else if (lrs == null) {
                    Thread.sleep(5000);
                } else if (lrs.size() == 0) {
                    lrs = writeQueue.take();
                    numPoisonPills++;
                } else {
                    lrs = writeQueue.take();

                    File cgf = makeTempGraphFile(cgout);

                    CortexGraphWriter cgw = makeGraphWriter(sampleName, kmerSize, cgf);

                    for (CortexRecord lr : lrs) {
                        cgw.addRecord(lr);
                    }

                    cgw.close();

                    log.info("  -  wrote: {} records to {} [{}]", lrs.size(), cgf.getAbsolutePath(), Thread.currentThread().getName());

                    graphs.add(new CortexGraph(cgf));
                }
            }
        } catch (InterruptedException e) {
            throw new CortexJDKException("Interrupted while sorting records", e);
        }

        return graphs;
    }
}
