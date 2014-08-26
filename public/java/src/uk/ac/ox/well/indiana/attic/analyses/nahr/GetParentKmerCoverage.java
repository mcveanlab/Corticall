package uk.ac.ox.well.indiana.attic.analyses.nahr;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.io.PrintStream;

public class GetParentKmerCoverage extends Module {
    @Argument(fullName="merged", shortName="m", doc="Merged sample, parents, and mask (in Cortex format)")
    public CortexGraph MERGED;

    @Output(fullName="oref0", shortName="o0", doc="Ref0 coverage")
    public PrintStream oref0;

    @Output(fullName="oref1", shortName="o1", doc="Ref1 coverage")
    public PrintStream oref1;

    @Override
    public void execute() {
        int ref0 = 0;
        int ref1 = 0;
        int refBoth = 0;
        int refTotal = 0;

        log.info("Computing coverage distribution for diagnostic kmers...");

        oref1.println("ref\tsample");
        oref0.println("ref\tsample");

        int numRecords = 0;
        for (CortexRecord cr : MERGED) {
            if (numRecords % (MERGED.getNumRecords() / 10) == 0) {
                log.info("  processed {}/{} (~{}%) records", numRecords, MERGED.getNumRecords(), String.format("%.2f", 100.0f*((float) numRecords / (float) MERGED.getNumRecords())));
            }
            numRecords++;

            int covSample = cr.getCoverage(0);
            int covRef0 = cr.getCoverage(1);
            int covRef1 = cr.getCoverage(2);
            int covMask = cr.getCoverage(3);

            if (covMask == 0) {
                if (covRef0 == 0 && covRef1 > 0) {
                    // HB3
                    ref1++;
                    oref1.println(covRef1 + "\t" + covSample);
                } else if (covRef0 > 0 && covRef1 == 0) {
                    // 3D7
                    ref0++;
                    oref0.println(covRef0 + "\t" + covSample);
                } else {
                    refBoth++;
                }
            }

            refTotal++;
        }

        log.info("Number of diagnostic kmers:");
        log.info("  ref0: {}/{} (~{}%)", ref0, refTotal, String.format("%.2f", 100.0f*((float) ref0)/((float) refTotal)));
        log.info("  ref1: {}/{} (~{}%)", ref1, refTotal, String.format("%.2f", 100.0f*((float) ref1)/((float) refTotal)));
        log.info("  refBoth: {}/{} (~{}%)", refBoth, refTotal, String.format("%.2f", 100.0f*((float) refBoth)/((float) refTotal)));
    }
}
