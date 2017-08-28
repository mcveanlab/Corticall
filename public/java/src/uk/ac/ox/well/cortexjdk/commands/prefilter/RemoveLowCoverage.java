package uk.ac.ox.well.cortexjdk.commands.prefilter;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraphWriter;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;

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

        for (CortexRecord cr : ROI) {
            if (cr.getCoverage(0) >= MIN_COVERAGE) {
                cgw.addRecord(cr);
            } else {
                lgw.addRecord(cr);
            }
        }

        cgw.close();
        lgw.close();
    }
}
