package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

public class KmerCoverageDistribution extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="color", shortName="c", doc="Color to use")
    public Integer COLOR = 0;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<Integer, Integer> kmerCoverageHist = new TreeMap<Integer, Integer>();

        log.info("Processing records...");
        int numRecords = 0;
        for (CortexRecord cr : GRAPH) {
            if (numRecords % (GRAPH.getNumRecords() / 10) == 0) {
                log.info("  {}/{} records", numRecords, GRAPH.getNumRecords());
            }

            int cov = cr.getCoverage(COLOR);

            if (!kmerCoverageHist.containsKey(cov)) {
                kmerCoverageHist.put(cov, 1);
            } else {
                kmerCoverageHist.put(cov, kmerCoverageHist.get(cov) + 1);
            }
        }

        out.println(Joiner.on("\t").join("cov", "count"));
        for (Integer cov : kmerCoverageHist.keySet()) {
            out.println(Joiner.on("\t").join(cov, kmerCoverageHist.get(cov)));
        }
    }
}
