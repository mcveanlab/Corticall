package uk.ac.ox.well.indiana.commands.gg;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraphWriter;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.math.MoreMathUtils;

import java.io.File;
import java.util.Map;

public class EmitNovelKmers extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="thresholds", shortName="t", doc="Thresholds")
    public File THRESHOLDS;

    @Output
    public File out;

    @Override
    public void execute() {
        log.info("Reading threshold information...");

        int p1Threshold = 0, p2Threshold = 0, bThreshold = 0;

        TableReader tr = new TableReader(THRESHOLDS);
        for (Map<String, String> te : tr) {
            String className = te.get("class");
            int threshold = Integer.valueOf(te.get("threshold"));

            if (className.equals("X10")) {
                p1Threshold = threshold;
            } else if (className.equals("X01")) {
                p2Threshold = threshold;
            } else if (className.equals("X11")) {
                bThreshold = threshold;
            }
        }

        int threshold = MoreMathUtils.max(p1Threshold, p2Threshold, bThreshold);

        log.info("  threshold: {} ({} {} {})", threshold, p1Threshold, p2Threshold, bThreshold);

        log.info("Processing graph...");
        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(GRAPH.getHeader());

        int numRecords = 0, numNovelRecords = 0;
        for (CortexRecord cr : GRAPH) {
            if (numRecords % (GRAPH.getNumRecords() / 10) == 0) {
                log.info("  {}/{} records ({} novel)", numRecords, GRAPH.getNumRecords(), numNovelRecords);
            }

            int cov    = cr.getCoverage(0);
            int cov_p1 = cr.getCoverage(1);
            int cov_p2 = cr.getCoverage(2);
            int cov_r1 = cr.getCoverage(3);
            int cov_r2 = cr.getCoverage(4);

            if (cov > threshold && cov_p1 == 0 && cov_p2 == 0 && cov_r1 == 0 && cov_r2 == 0) {
                cgw.addRecord(cr);
                numNovelRecords++;
            }

            numRecords++;
        }

        cgw.close();

        log.info("  {} total novel records", numNovelRecords);
    }
}
