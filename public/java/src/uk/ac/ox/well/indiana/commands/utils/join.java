package uk.ac.ox.well.indiana.commands.utils;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.collection.CortexCollection;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraphWriter;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.math.MoreMathUtils;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;

import java.io.File;

public class join extends Module {
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
    }
}
