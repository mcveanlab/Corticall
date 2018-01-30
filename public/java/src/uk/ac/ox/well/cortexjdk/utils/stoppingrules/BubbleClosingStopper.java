package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalState;

/**
 * Created by kiran on 08/05/2017.
 */
public class BubbleClosingStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    @Override
    public boolean hasTraversalSucceeded(TraversalState<CortexVertex> s) {
        //return s.getPreviousGraph().containsVertex(s.getCurrentVertex());

        return false;
    }

    @Override
    public boolean hasTraversalFailed(TraversalState<CortexVertex> s) {
        return s.getCurrentGraphSize() > 10000 || s.getCurrentJunctionDepth() >= 2 || s.getNumAdjacentEdges() == 0;
    }
}
