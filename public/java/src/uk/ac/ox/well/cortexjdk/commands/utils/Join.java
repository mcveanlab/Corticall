package uk.ac.ox.well.cortexjdk.commands.utils;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexCollection;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraphWriter;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;

import java.io.File;
import java.util.ArrayList;

public class Join extends Module {
    @Argument(fullName="graphs", shortName="g", doc="Graphs")
    public ArrayList<CortexGraph> GRAPHS;

    @Output
    public File out;

    @Override
    public void execute() {
        CortexCollection cc = new CortexCollection(GRAPHS);

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(cc.getHeader());

        log.info("Joining graphs:");
        for (int c = 0; c < cc.getNumColors(); c++) {
            log.info("  {}: {} ({} kmers)", c, cc.getSampleName(c), cc.getGraph(c).getNumRecords());
        }

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing graph...")
                .message("records processed")
                .updateRecord(cc.getGraph(0).getNumRecords() / 10)
                .make(log);

        for (CortexRecord cr : cc) {
            cgw.addRecord(cr);

            pm.update();
        }

        cgw.close();
    }
}
