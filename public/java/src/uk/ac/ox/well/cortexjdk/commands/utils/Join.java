package uk.ac.ox.well.cortexjdk.commands.utils;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.collection.CortexCollection;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraphWriter;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;

import java.io.File;

public class Join extends Module {
    @Argument(fullName="graphs", shortName="g", doc="Graphs")
    public CortexCollection GRAPHS;

    @Output
    public File out;

    @Override
    public void execute() {
        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(GRAPHS.getHeader());

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing graph...")
                .message("records processed (estimated)")
                .updateRecord(GRAPHS.getGraph(0).getNumRecords() / 10)
                .maxRecord(GRAPHS.getGraph(0).getNumRecords())
                .make(log);

        for (CortexRecord cr : GRAPHS) {
            cgw.addRecord(cr);

            pm.update();
        }

        cgw.close();
    }
}
