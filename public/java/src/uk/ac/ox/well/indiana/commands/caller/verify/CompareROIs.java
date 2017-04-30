package uk.ac.ox.well.indiana.commands.caller.verify;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

public class CompareROIs extends Module {
    @Argument(fullName="expected", shortName="e", doc="Expected")
    public CortexGraph EXPECTED;

    @Argument(fullName="actual", shortName="a", doc="Actual")
    public CortexGraph ACTUAL;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CortexKmer, CortexRecord> ars = new HashMap<>();
        Map<CortexKmer, Boolean> arsUsed = new HashMap<>();

        for (CortexRecord ar : ACTUAL) {
            ars.put(ar.getCortexKmer(), ar);
            arsUsed.put(ar.getCortexKmer(), false);
        }

        int privateToExpected = 0, privateToActual = 0, overlap = 0;

        for (CortexRecord er : EXPECTED) {
            boolean presentInActual = ars.containsKey(er.getCortexKmer());
            CortexRecord ar = ars.get(er.getCortexKmer());

            log.info("exp={} act={} er={} ar={}", true, presentInActual, er, ar);
            out.println("exp=" + true + " act=" + presentInActual + " er=" + er + " ar=" + ar);

            arsUsed.put(er.getCortexKmer(), true);

            if (presentInActual) {
                overlap++;
            } else {
                privateToExpected++;
            }
        }

        for (CortexKmer ak : ars.keySet()) {
            if (!arsUsed.containsKey(ak)) {
                CortexRecord ar = ars.get(ak);

                log.info("exp={} act={} er={} ar={}", false, true, null, ar);
                out.println("exp=" + false + " act=" + true + " er=null ar=" + ar);
            }
        }

        log.info("pe={} pa={} o={}", privateToExpected, privateToActual, overlap);
        out.println("pe=" + privateToExpected + " pa=" + privateToActual + " o=" + overlap);
    }
}
