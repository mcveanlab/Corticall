package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataTable;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

public class ConditionalKmerCoverageDistribution extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Output
    public PrintStream out;

    private void incrementElement(Map<Integer, Map<String, Integer>> covCounts, int cov, String code) {
        if (!covCounts.containsKey(cov)) {
            covCounts.put(cov, new HashMap<String, Integer>());
        }

        if (!covCounts.get(cov).containsKey(code)) {
            covCounts.get(cov).put(code, 0);
        }

        covCounts.get(cov).put(code, covCounts.get(cov).get(code) + 1);
    }

    private int getElement(Map<Integer, Map<String, Integer>> covCounts, int cov, String code) {
        if (covCounts.containsKey(cov) && covCounts.get(cov).containsKey(code)) {
            return covCounts.get(cov).get(code);
        }

        return 0;
    }

    @Override
    public void execute() {
        Map<Integer, Map<String, Integer>> covCounts = new TreeMap<Integer, Map<String, Integer>>();

        for (CortexRecord cr : GRAPH) {
            int cov    = cr.getCoverage(0);
            int cov_p1 = cr.getCoverage(1);
            int cov_p2 = cr.getCoverage(2);
            int cov_r1 = cr.getCoverage(3);
            int cov_r2 = cr.getCoverage(4);

            if        (cov_p1 >  0 && cov_p2 == 0 && cov_r1 == 0 && cov_r2 == 0) {
                incrementElement(covCounts, cov, "1000");
            } else if (cov_p1 == 0 && cov_p2 >  0 && cov_r1 == 0 && cov_r2 == 0) {
                incrementElement(covCounts, cov, "0100");
            } else if (cov_p1 == 0 && cov_p2 == 0 && cov_r1 >  0 && cov_r2 == 0) {
                incrementElement(covCounts, cov, "0010");
            } else if (cov_p1 == 0 && cov_p2 == 0 && cov_r1 == 0 && cov_r2 >  0) {
                incrementElement(covCounts, cov, "0001");
            }
        }

        out.println(Joiner.on("\t").join("", "1000", "0100", "0010", "0001"));

        for (Integer cov : covCounts.keySet()) {
            int count_p1 = getElement(covCounts, cov, "1000");
            int count_p2 = getElement(covCounts, cov, "0100");
            int count_r1 = getElement(covCounts, cov, "0010");
            int count_r2 = getElement(covCounts, cov, "0001");

            out.println(Joiner.on("\t").join(cov, count_p1, count_p2, count_r1, count_r2));
        }
    }
}
