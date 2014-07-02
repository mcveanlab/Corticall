package uk.ac.ox.well.indiana.attic.analyses.nahr;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexMap;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;

import java.io.PrintStream;

public class GetParentKmerCoverage extends Module {
    @Argument(fullName="parents", shortName="p", doc="Parents (in Cortex format)")
    public CortexGraph PARENTS;

    @Argument(fullName="sample", shortName="s", doc="Sample (in Cortex format)", required=false)
    public CortexMap SAMPLE;

    @Argument(fullName="maskedKmers", shortName="m", doc="Masked kmers", required=false)
    public CortexMap MASKED_KMERS;

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

        if (SAMPLE != null) {
            oref1.println("ref\tsample");
            oref0.println("ref\tsample");
        }

        int numRecords = 0;
        for (CortexRecord cr : PARENTS) {
            if (numRecords % (PARENTS.getNumRecords() / 5) == 0) {
                log.info("  processed {}/{} (~{}%) records", numRecords, PARENTS.getNumRecords(), String.format("%.2f", 100.0f*((float) numRecords / (float) PARENTS.getNumRecords())));
            }
            numRecords++;

            if (MASKED_KMERS == null || !MASKED_KMERS.containsKey(cr.getKmer())) {
                if (cr.getCoverage(0) == 0 && cr.getCoverage(1) > 0) {
                    // HB3
                    ref1++;

                    if (SAMPLE != null) {
                        int cov = SAMPLE.containsKey(cr.getKmer()) ? SAMPLE.get(cr.getKmer()).getCoverage(0) : 0;
                        oref1.println(cr.getCoverage(1) + "\t" + cov);
                    } else {
                        oref1.println(cr.getCoverage(1));
                    }
                } else if (cr.getCoverage(0) > 0 && cr.getCoverage(1) == 0) {
                    // 3D7
                    ref0++;

                    if (SAMPLE != null) {
                        int cov = SAMPLE.containsKey(cr.getKmer()) ? SAMPLE.get(cr.getKmer()).getCoverage(0) : 0;
                        oref0.println(cr.getCoverage(0) + "\t" + cov);
                    } else {
                        oref0.println(cr.getCoverage(0));
                    }
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
