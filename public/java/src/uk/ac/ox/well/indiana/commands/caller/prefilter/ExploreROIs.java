package uk.ac.ox.well.indiana.commands.caller.prefilter;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngineFactory;
import uk.ac.ox.well.indiana.utils.visualizer.GraphVisualizer;

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
        GraphVisualizer gv = new GraphVisualizer(9000);

        for (CortexRecord rr : ROIS) {
            CortexRecord cr = GRAPH.findRecord(rr.getCortexKmer());

            log.info("{} {} {} {} {} {}",
                    ROIS.getSampleName(8),
                    ROIS.getSampleName(9),
                    ROIS.getSampleName(10),
                    ROIS.getSampleName(11),
                    rr,
                    cr
            );

            //TraversalEngine e = new TraversalEngineFactory()
            //DirectedGraph<CortexVertex, CortexEdge> g =
        }
    }
}
