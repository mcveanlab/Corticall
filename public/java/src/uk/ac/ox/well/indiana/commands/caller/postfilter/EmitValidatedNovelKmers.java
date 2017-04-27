package uk.ac.ox.well.indiana.commands.caller.postfilter;

import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.*;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;

import java.io.File;

/**
 * Created by kiran on 26/04/2017.
 */
public class EmitValidatedNovelKmers extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Output
    public File out;

    @Override
    public void execute() {
        /*
        Colour 0: sample name: '803' (drafts)
        Colour 1: sample name: 'PG0443-C' (clean)
        Colour 2: sample name: 'PG0443-C' (tips)
        Colour 3: sample name: 'PG0443-C' (raw)
        Colour 4: sample name: 'GB4' (drafts)
        Colour 5: sample name: 'PG0050-CX2' (clean)
        Colour 6: sample name: 'PG0050-CX2' (tips)
        Colour 7: sample name: 'PG0050-CX2' (raw)
        Colour 8: sample name: '34F5' (drafts)
        Colour 9: sample name: 'PG0476-C' (clean)
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
        //cgw.setHeader(makeCortexHeader());

        for (CortexRecord cr : GRAPH) {
            if (cr.getCoverage(0) == 0 && cr.getCoverage(1) == 0 && cr.getCoverage(2) == 0 && cr.getCoverage(3) == 0 &&
                cr.getCoverage(4) == 0 && cr.getCoverage(5) == 0 && cr.getCoverage(6) == 0 && cr.getCoverage(7) == 0 &&
                cr.getCoverage(8)  > 0 && cr.getCoverage(9)  > 0 && cr.getCoverage(10) > 0 && cr.getCoverage(11) > 0) {

                cgw.addRecord(cr);
            }

            pm.update();
        }
    }

    @NotNull
    private CortexHeader makeCortexHeader() {
        CortexHeader ch = new CortexHeader();
        ch.setVersion(6);
        ch.setNumColors(1);
        ch.setKmerSize(GRAPH.getKmerSize());
        ch.setKmerBits(GRAPH.getKmerBits());

        CortexColor cc = new CortexColor();
        cc.setCleanedAgainstGraph(false);
        cc.setCleanedAgainstGraphName("");
        cc.setErrorRate(0);
        cc.setLowCovgKmersRemoved(false);
        cc.setLowCovgSupernodesRemoved(false);
        cc.setTipClippingApplied(false);
        cc.setTotalSequence(0);
        cc.setSampleName(GRAPH.getSampleName(8));

        ch.addColor(cc);

        return ch;
    }
}
