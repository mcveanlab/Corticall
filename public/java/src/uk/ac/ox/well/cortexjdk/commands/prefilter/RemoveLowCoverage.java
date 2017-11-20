package uk.ac.ox.well.cortexjdk.commands.prefilter;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraphWriter;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;

import java.io.File;

/**
 * Created by kiran on 28/08/2017.
 */
public class RemoveLowCoverage extends Module {
    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="minCoverage", shortName="m", doc="Coverage limit")
    public Integer MIN_COVERAGE = 10;

    @Output
    public File out;

    @Output(fullName="lowcoverage_out", shortName="co", doc="Low coverage output file")
    public File lowcoverage_out;

    @Override
    public void execute() {
        CortexGraphWriter cgw = new CortexGraphWriter(out);
        CortexGraphWriter lgw = new CortexGraphWriter(lowcoverage_out);

        cgw.setHeader(ROI.getHeader());
        lgw.setHeader(ROI.getHeader());

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Removing low-coverage records (< " + MIN_COVERAGE + ")")
                .message("records processed")
                .maxRecord(ROI.getNumRecords())
                .make(log);

        int numKept = 0, numExcluded = 0;
        for (CortexRecord cr : ROI) {
            if (cr.getCoverage(0) >= MIN_COVERAGE) {
                cgw.addRecord(cr);
                numKept++;
            } else {
                lgw.addRecord(cr);
                numExcluded++;
            }

            pm.update();
        }

        cgw.close();
        lgw.close();

        log.info("  {}/{} ({}%) kept, {}/{} ({}%) excluded",
                numKept,     ROI.getNumRecords(), 100.0f * (float) numKept / (float) ROI.getNumRecords(),
                numExcluded, ROI.getNumRecords(), 100.0f * (float) numExcluded / (float) ROI.getNumRecords()
        );
    }
}
