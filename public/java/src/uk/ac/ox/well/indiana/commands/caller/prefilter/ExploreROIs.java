package uk.ac.ox.well.indiana.commands.caller.prefilter;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.stoppingconditions.ExplorationStopper;
import uk.ac.ox.well.indiana.utils.traversal.*;
import uk.ac.ox.well.indiana.utils.visualizer.GraphVisualizer;

import java.util.*;

/**
 * Created by kiran on 30/05/2017.
 */
public class ExploreROIs extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="rois", shortName="r", doc="ROIs")
    public CortexGraph ROIS;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(ROIS.getSampleName(11));
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(Arrays.asList("PG0443-C", "PG0050-CX2"));

        log.info("Colors:");
        log.info(" -   child: {}", childColor);
        log.info(" - parents: {}", parentColors);

        Set<Integer> displayColors = new TreeSet<>();
        displayColors.add(childColor);
        displayColors.addAll(parentColors);

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(childColor)
                .joiningColors(parentColors)
                .displayColors(displayColors)
                .traversalDirection(TraversalEngineConfiguration.TraversalDirection.BOTH)
                .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.OR)
                .stopper(ExplorationStopper.class)
                .graph(GRAPH)
                .rois(ROIS)
                .make();

        GraphVisualizer gv = new GraphVisualizer(9000);

        List<CortexRecord> rrs = new ArrayList<>();
        for (CortexRecord rr : ROIS) {
            rrs.add(rr);
        }

        for (CortexRecord rr : rrs) {
            CortexRecord cr = GRAPH.findRecord(rr.getCortexKmer());

            log.info("{} {} {} {} {} {}",
                    ROIS.getSampleName(8),
                    ROIS.getSampleName(9),
                    ROIS.getSampleName(10),
                    ROIS.getSampleName(11),
                    rr,
                    cr
            );

            DirectedGraph<CortexVertex, CortexEdge> g = e.dfs(rr.getKmerAsString());

            gv.display(g, displayColors);

            log.info("Press enter for next graph");
            System.console().readLine();
        }
    }
}
