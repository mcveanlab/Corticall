package uk.ac.ox.well.indiana.commands.utils;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraphWriter;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.util.Arrays;

public class sort extends Module {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Output
    public CortexGraphWriter out;

    @Override
    public void execute() {
        CortexRecord[] records = new CortexRecord[(int) CORTEX_GRAPH.getNumRecords()];

        log.info("Processing Cortex records...");
        int recordsProcessed = 0;
        for (CortexRecord cr : CORTEX_GRAPH) {
            records[recordsProcessed] = cr;

            recordsProcessed++;
            if (recordsProcessed % (CORTEX_GRAPH.getNumRecords() / 10) == 0) {
                log.info("  {}/{} records", recordsProcessed, CORTEX_GRAPH.getNumRecords());
            }
        }

        log.info("Sorting records (this may take a while)...");
        Arrays.sort(records, 0, records.length);

        log.info("Writing records...");
        recordsProcessed = 0;
        out.setHeader(CORTEX_GRAPH.getHeader());
        for (CortexRecord cr : records) {
            recordsProcessed++;

            out.addRecord(cr);

            if (recordsProcessed % (records.length / 10) == 0) {
                log.info("  {}/{} records", recordsProcessed, records.length);
            }
        }
    }
}
