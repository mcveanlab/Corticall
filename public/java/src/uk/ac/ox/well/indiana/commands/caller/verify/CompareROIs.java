package uk.ac.ox.well.indiana.commands.caller.verify;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class CompareROIs extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="truth", shortName="t", doc="Truth")
    public CortexGraph TRUTH;

    @Argument(fullName="eval", shortName="e", doc="Eval")
    public CortexGraph EVAL;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CortexKmer, CortexRecord> trs = new HashMap<>();
        Map<CortexKmer, CortexRecord> ers = new HashMap<>();
        Set<CortexKmer> all = new HashSet<>();

        for (CortexRecord tr : TRUTH) {
            trs.put(tr.getCortexKmer(), tr);
            all.add(tr.getCortexKmer());
        }

        for (CortexRecord er : EVAL) {
            ers.put(er.getCortexKmer(), er);
            all.add(er.getCortexKmer());
        }

        int privateToTruth = 0, privateToEval = 0, overlap = 0;

        for (CortexKmer ck : all) {
            if (trs.containsKey(ck)) {
                if (ers.containsKey(ck)) {
                    overlap++;
                } else {
                    privateToTruth++;

                    log.debug("{} {}", TRUTH.getSampleName(0), GRAPH.findRecord(ck));
                }
            } else {
                privateToEval++;
            }
        }

        log.info("t={} e={} pt={} pe={} o={}", trs.size(), ers.size(), privateToTruth, privateToEval, overlap);
        out.println("t=" + trs.size() + " e=" + ers.size() + " pt=" + privateToTruth + " pe=" + privateToEval + " o=" + overlap);
    }
}
