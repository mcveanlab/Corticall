package uk.ac.ox.well.cortexjdk.playground.eval;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;

import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;

/**
 * Created by kiran on 27/08/2017.
 */
public class SummarizeFilteredNovelKmers extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public ArrayList<CortexGraph> GRAPHS;

    @Override
    public void execute() {
        Map<Long, String> counts = new TreeMap<>();

        for (CortexGraph cg : GRAPHS) {
            String filterName = cg.getCortexFile().getName();
            filterName = filterName.replaceFirst(".+rois.", "");
            filterName = filterName.replaceAll("ctx", "");

            counts.put(cg.getNumRecords(), filterName);
        }

        for (long l : counts.keySet()) {
            String[] s = counts.get(l).split("\\.");

            if (s.length > 1) {
                log.info("{} {}", s[0], l);
            } else {
                log.info("{} {}", "unfiltered", l);
            }
        }
    }
}
