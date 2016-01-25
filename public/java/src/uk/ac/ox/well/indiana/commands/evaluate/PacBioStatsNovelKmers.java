package uk.ac.ox.well.indiana.commands.evaluate;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.io.PrintStream;

public class PacBioStatsNovelKmers extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int tp = 0, tn = 0, fp = 0, fn = 0, un = 0;

        log.info("Processing records...");
        int index = 0;
        for (CortexRecord cr : GRAPH) {
            if (index % (GRAPH.getNumRecords() / 10) == 0) {
                log.info("  {}/{} records", index, GRAPH.getNumRecords());
            }
            index++;

            boolean calledNovel = cr.getCoverage(0) > 0;
            boolean rejectedNovel = cr.getCoverage(1) > 0;
            boolean presentInPacBio = cr.getCoverage(2) > 0;
            boolean presentInChild = cr.getCoverage(3) > 0;
            boolean presentInMom = cr.getCoverage(4) > 0;
            boolean presentInDad = cr.getCoverage(5) > 0;
            boolean presentInParent = presentInMom || presentInDad;

            // Called a novel
            if (calledNovel) {
                // Correct? TP
                if (presentInPacBio && presentInChild && !presentInParent) {
                    tp++;

                    log.debug("tp: {}", cr);
                } else {
                    fp++;

                    log.debug("fp: {}", cr);
                }
            } else if (rejectedNovel) {
                // Correct? TN
                if (!presentInPacBio) {
                    tn++;

                    log.debug("tn: {}", cr);
                } else {
                    fn++;

                    log.debug("fn: {}", cr);
                }
            } else {
                if (presentInPacBio && presentInChild && !presentInParent) {
                    fn++;

                    log.info("huh?: {}", cr);
                }
            }

        }

        float sens = (float) tp / (float) (tp + fn);
        float spec = (float) tn / (float) (tn + fp);

        out.println("tp: " + tp);
        out.println("tn: " + tn);
        out.println("fp: " + fp);
        out.println("fn: " + fn);
        out.println("un: " + un);
        out.println("sens: " + sens);
        out.println("spec: " + spec);
    }
}
