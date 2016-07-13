package uk.ac.ox.well.indiana.attic.gg;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraphWriter;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;

import java.io.File;

public class EmitNovelKmers extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="threshold", shortName="t", doc="Threshold file", required=false)
    public File THRESHOLD_FILE;

    @Argument(fullName="upperThreshold", shortName="u", doc="Upper threshold")
    public Integer UPPER_THRESHOLD = 1000000;

    @Output
    public File out;

    @Override
    public void execute() {
        int threshold = 0;

        if (THRESHOLD_FILE != null) {
            LineReader lr = new LineReader(THRESHOLD_FILE);
            threshold = Integer.valueOf(lr.getNextRecord().trim());
        }

        log.info("Processing graph...");
        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(GRAPH.getHeader());

        int numRecords = 0, numNovelRecords = 0;
        for (CortexRecord cr : GRAPH) {
            if (GRAPH.getNumRecords() > 10 && numRecords % (GRAPH.getNumRecords() / 10) == 0) {
                log.info("  {}/{} records ({} novel)", numRecords, GRAPH.getNumRecords(), numNovelRecords);
            }

            int cov = cr.getCoverage(0);

            if (cr.getCortexKmer().equals(new CortexKmer("AAATAAATAAAAAATATAAATATAAATATAAATATATATATATATAT"))) {
                log.info("Hi!");
            }

            boolean hasCoverageInOtherColors = false;
            for (int c = 1; c < cr.getNumColors(); c++) {
                if (cr.getCoverage(c) > 0) {
                    hasCoverageInOtherColors = true;
                    break;
                }
            }

            if (cov > threshold && cov < UPPER_THRESHOLD && !hasCoverageInOtherColors) {
                cgw.addRecord(cr);
                numNovelRecords++;
            }

            numRecords++;
        }

        cgw.close();

        log.info("  {} total novel records", numNovelRecords);
    }
}
