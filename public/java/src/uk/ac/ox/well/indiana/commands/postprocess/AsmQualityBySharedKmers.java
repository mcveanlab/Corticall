package uk.ac.ox.well.indiana.commands.postprocess;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class AsmQualityBySharedKmers extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="eval", shortName="e", doc="Evaluation sample")
    public String EVAL;

    @Argument(fullName="val", shortName="v", doc="Validation sample(s)")
    public ArrayList<String> VALS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int evalColor = GRAPH.getColorForSampleName(EVAL);
        List<Integer> valColors = new ArrayList<>();
        for (String val : VALS) {
            valColors.add(GRAPH.getColorForSampleName(val));
        }

        int numKmers = 0;
        int numWrongKmers = 0;

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Examining records")
                .message("records examined")
                .maxRecord(GRAPH.getNumRecords())
                .make(log);

        for (CortexRecord cr : GRAPH) {
            boolean isInEval = cr.getCoverage(evalColor) > 0;
            boolean isInVal = false;

            for (int c : valColors) {
                if (cr.getCoverage(c) > 0) {
                    isInVal = true;
                    break;
                }
            }

            if (isInEval && !isInVal) {
                numWrongKmers++;
            }
            numKmers++;

            pm.update();
        }

        log.info("numKmers={} numWrongKmers={} Q={}", numKmers, numWrongKmers, -10.0*Math.log10((double) numWrongKmers / (double) numKmers));
    }
}
