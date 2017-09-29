package uk.ac.ox.well.cortexjdk.commands.call.verify;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraphWriter;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;

import java.io.File;

/**
 * Created by kiran on 26/04/2017.
 */
public class EmitValidatedNovelKmers extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="covLimit", shortName="c", doc="Coverage limit")
    public Integer COV_LIMIT = 0;

    @Output
    public File out;

    @Override
    public void execute() {
        /*
        Colour 0: sample name: '803' (drafts)
        Colour 1: sample name: 'PG0443-C' (graph)
        Colour 2: sample name: 'PG0443-C' (tips)
        Colour 3: sample name: 'PG0443-C' (raw)
        Colour 4: sample name: 'GB4' (drafts)
        Colour 5: sample name: 'PG0050-CX2' (graph)
        Colour 6: sample name: 'PG0050-CX2' (tips)
        Colour 7: sample name: 'PG0050-CX2' (raw)
        Colour 8: sample name: '34F5' (drafts)
        Colour 9: sample name: 'PG0476-C' (graph)
        Colour 10: sample name: 'PG0476-C' (tips)
        Colour 11: sample name: 'PG0476-C' (raw)
        */

        ProgressMeter pm = new ProgressMeterFactory()
                .maxRecord(GRAPH.getNumRecords())
                .header("Processing records")
                .message("records processed")
                .make(log);


        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(GRAPH.getHeader());

        for (CortexRecord cr : GRAPH) {
            if (cr.getCoverage(0) == 0 && cr.getCoverage(1) == 0 && cr.getCoverage(2) == 0 && cr.getCoverage(3) == 0 &&
                cr.getCoverage(4) == 0 && cr.getCoverage(5) == 0 && cr.getCoverage(6) == 0 && cr.getCoverage(7) == 0 &&
                cr.getCoverage(8)  > 0 && cr.getCoverage(9)  > 0 && cr.getCoverage(10) > 0 && cr.getCoverage(11) > COV_LIMIT) {

                cgw.addRecord(cr);
            }

            pm.update();
        }

        cgw.close();
    }
}
