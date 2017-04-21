package uk.ac.ox.well.indiana.commands.postprocess;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
        int numMissedKmers = 0;

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Examining records")
                .message("records examined")
                .maxRecord(GRAPH.getNumRecords())
                .make(log);

        Set<CortexKmer> wrongKmers = new HashSet<>();

        for (CortexRecord cr : GRAPH) {
            boolean isInEval = cr.getCoverage(evalColor) > 0;
            boolean isInVal = false;

            for (int c : valColors) {
                if (cr.getCoverage(c) > 0) {
                    isInVal = true;
                    break;
                }
            }

            if (isInEval || isInVal) {
                numKmers++;

                if (isInEval && !isInVal) {
                    numWrongKmers++;

                    wrongKmers.add(cr.getCortexKmer());
                } else if (!isInEval && isInVal) {
                    numMissedKmers++;
                }
            }

            pm.update();
        }

        /*
        int numEvents = 0;
        Set<CortexKmer> usedKmers = new HashSet<>();
        for (CortexKmer ck : wrongKmers) {
            if (!usedKmers.contains(ck)) {
                usedKmers.add(ck);

                String stretch = CortexUtils.getSeededStretch(GRAPH, ck.getKmerAsString(), evalColor, false);
                for (int i = 0; i <= stretch.length() - GRAPH.getKmerSize(); i++) {
                    CortexKmer nk = new CortexKmer(stretch.substring(i, i + GRAPH.getKmerSize()));

                    usedKmers.add(nk);
                }

                numEvents++;
            }
        }
        */

        log.info("numKmers={} numWrongKmers={} numMissedKmers={} Q={} Q(/10)={} Q(/k)={}",
                numKmers,
                numWrongKmers,
                numMissedKmers,
                -10.0*Math.log10((double) numWrongKmers / (double) numKmers),
                -10.0*Math.log10((((double) numWrongKmers)/10.0) / (double) numKmers),
                -10.0*Math.log10((((double) numWrongKmers)/ (double) GRAPH.getKmerSize()) / (double) numKmers)
        );
    }
}
