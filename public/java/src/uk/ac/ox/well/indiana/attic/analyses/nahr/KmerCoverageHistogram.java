package uk.ac.ox.well.indiana.attic.analyses.nahr;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

public class KmerCoverageHistogram extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public HashMap<String, CortexGraph> GRAPHS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<Integer, Map<String, Integer>> covHist = new TreeMap<Integer, Map<String, Integer>>();

        for (String label : GRAPHS.keySet()) {
            CortexGraph cg = GRAPHS.get(label);

            for (CortexRecord cr : cg) {
                int cov = cr.getCoverage(0);

                if (!covHist.containsKey(cov)) {
                    covHist.put(cov, new TreeMap<String, Integer>());
                }

                if (!covHist.get(cov).containsKey(label)) {
                    covHist.get(cov).put(label, 1);
                } else {
                    covHist.get(cov).put(label, 1 + covHist.get(cov).get(label));
                }
            }
        }

        TableWriter tw = new TableWriter(out);
        for (int cov : covHist.keySet()) {
            Map<String, String> te = new HashMap<String, String>();
            te.put("cov", String.valueOf(cov));

            for (String label : GRAPHS.keySet()) {
                int count = (covHist.get(cov).containsKey(label)) ? covHist.get(cov).get(label) : 0;

                te.put(label, String.valueOf(count));
            }

            tw.addEntry(te);
        }
    }
}
