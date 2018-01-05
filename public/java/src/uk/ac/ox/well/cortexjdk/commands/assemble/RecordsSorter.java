package uk.ac.ox.well.cortexjdk.commands.assemble;

import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import org.slf4j.Logger;

public class RecordsSorter implements Runnable {
    private BlockingQueue<List<CortexRecord>> readQueue;
    private BlockingQueue<List<CortexRecord>> writeQueue;
    private Logger log;

    public RecordsSorter(BlockingQueue<List<CortexRecord>> readQueue, BlockingQueue<List<CortexRecord>> writeQueue, Logger log) {
        this.readQueue = readQueue;
        this.writeQueue = writeQueue;
        this.log = log;
    }

    @Override
    public void run() {
        try {
            List<CortexRecord> lrs;

            while (true) {
                lrs = readQueue.peek();

                if (lrs == null) {
                    Thread.sleep(5000);
                } else if (lrs.size() == 0) {
                    break;
                } else {
                    lrs = readQueue.take();

                    Collections.sort(lrs);

                    log.info("  sorted {} records [{}]", lrs.size(), Thread.currentThread().getName());

                    writeQueue.put(lrs);
                }
            }

            writeQueue.put(new ArrayList<>());
        } catch (InterruptedException e) {
            throw new CortexJDKException("Interrupted while sorting records", e);
        }
    }
}
