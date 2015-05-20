package uk.ac.ox.well.indiana.commands.gg;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraphWriter;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;

import java.io.File;

public class EmitConfidentNovelKmers extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="threshold", shortName="t", doc="Threshold file")
    public File THRESHOLD_FILE;

    @Output
    public File out;

    @Override
    public void execute() {
        LineReader lr = new LineReader(THRESHOLD_FILE);
        int threshold = Integer.valueOf(lr.getNextRecord());

        log.info("Processing graph...");
        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(GRAPH.getHeader());

        int numRecords = 0, numNovelRecords = 0;
        for (CortexRecord cr : GRAPH) {
            if (numRecords % (GRAPH.getNumRecords() / 10) == 0) {
                log.info("  {}/{} records ({} novel)", numRecords, GRAPH.getNumRecords(), numNovelRecords);
            }

            int cov = cr.getCoverage(0);

            boolean hasCoverageInOtherColors = false;
            for (int c = 1; c < cr.getNumColors(); c++) {
                if (cr.getCoverage(c) > 0) {
                    hasCoverageInOtherColors = true;
                    break;
                }
            }

            if (cov > threshold && !hasCoverageInOtherColors) {
                cgw.addRecord(cr);
                numNovelRecords++;
            }

            numRecords++;
        }

        cgw.close();

        log.info("  {} total novel records", numNovelRecords);
    }
}
