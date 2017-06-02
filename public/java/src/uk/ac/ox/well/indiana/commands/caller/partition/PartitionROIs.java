package uk.ac.ox.well.indiana.commands.caller.partition;

import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DirectedWeightedMultigraph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.stoppingconditions.PartitionStopper;
import uk.ac.ox.well.indiana.utils.traversal.*;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.AND;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

public class PartitionROIs extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TraversalEngine e = new TraversalEngineFactory()
                .combinationOperator(OR)
                .traversalDirection(BOTH)
                .traversalColor(GRAPH.getColorForSampleName(CHILD))
                .joiningColors(GRAPH.getColorsForSampleNames(PARENTS))
                .recruitmentColors(GRAPH.getColorsForSampleNames(PARENTS))
                .rois(ROI)
                .stopper(PartitionStopper.class)
                .graph(GRAPH)
                .make();

        Set<CortexRecord> rois = new HashSet<>();
        Map<CortexKmer, Integer> p = new HashMap<>();
        Map<Integer, DirectedGraph<CortexVertex, CortexEdge>> m = new TreeMap<>();

        for (CortexRecord cr : ROI) {
            rois.add(cr);
            p.put(cr.getCortexKmer(), -1);
        }

        int group = 0;
        for (CortexKmer ck : p.keySet()) {
            log.info("{}", ck);

            if (p.containsKey(ck) && p.get(ck) == -1) {
                DirectedWeightedMultigraph<CortexVertex, CortexEdge> dfs = e.dfs(ck.getKmerAsString());

                if (dfs != null) {
                    int numNovels = 0;
                    for (CortexVertex cv : dfs.vertexSet()) {
                        if (p.containsKey(cv.getCr().getCortexKmer()) && p.get(cv.getCr().getCortexKmer()) == -1) {
                            p.put(cv.getCr().getCortexKmer(), group);
                            numNovels++;
                        }
                    }

                    log.info("  - {} {} {}", group, dfs.vertexSet().size(), numNovels);

                    m.put(group, dfs);

                    String contig = e.getContig(dfs, ck.getKmerAsString(), GRAPH.getColorForSampleName(CHILD));

                    out.println(">" + group);
                    out.println(contig);

                    group++;
                }
            }
        }

        log.info("Groups: {}", group);
    }
}
