package uk.ac.ox.well.cortexjdk.commands.quality;

import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.collection.CortexCollection;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.stoppingconditions.BubbleClosingStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingconditions.BubbleOpeningStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.AND;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;

/**
 * Created by kiran on 19/08/2017.
 */
public class AssessAssemblyQuality extends Module {
    @Argument(fullName="eval", shortName="e", doc="Graph to evaluate")
    public CortexGraph EVAL;

    @Argument(fullName="comp", shortName="c", doc="Graph to compare against (i.e. truth)")
    public CortexGraph COMP;

    @Override
    public void execute() {
        CortexCollection cc = new CortexCollection(COMP, EVAL);

        TraversalEngine eOpen = new TraversalEngineFactory()
                .traversalColor(0)
                .joiningColors(1)
                .traversalDirection(BOTH)
                .combinationOperator(AND)
                .stoppingRule(BubbleOpeningStopper.class)
                .graph(cc)
                .make();

        TraversalEngine eClose = new TraversalEngineFactory()
                .traversalColor(1)
                .joiningColors(0)
                .traversalDirection(BOTH)
                .combinationOperator(OR)
                .stoppingRule(BubbleClosingStopper.class)
                .graph(cc)
                .make();

        for (CortexRecord cr : cc) {
            if (isExclusiveToComp(cr) && isSinglyConnected(cr, 0)) {
                log.info("{}", cr);

                DirectedWeightedPseudograph<CortexVertex, CortexEdge> gOpen = eOpen.dfs(cr.getKmerAsString());

                for (CortexVertex cv : gOpen.vertexSet()) {
                    if (gOpen.inDegreeOf(cv) == 0 || gOpen.outDegreeOf(cv) == 0) {
                        eClose.getConfiguration().setPreviousTraversal(gOpen);

                        DirectedWeightedPseudograph<CortexVertex, CortexEdge> gClose = eClose.dfs(cv.getSk());
                    }
                }

            }
        }
    }

    private boolean isExclusiveToComp(CortexRecord cr) {
        return cr.getCoverage(0) > 0 && cr.getCoverage(1) == 0;
    }

    private boolean isSinglyConnected(CortexRecord cr, int c) {
        return cr.getInDegree(c) == 1 && cr.getOutDegree(c) == 1;
    }
}
