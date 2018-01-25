package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalState;

/**
 * Created by kiran on 14/05/2017.
 */
public class CycleCollapsingContigStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    @Override
    public boolean hasTraversalSucceeded(TraversalState<CortexVertex> s) {
        return s.getNumAdjacentEdges() == 0;
    }

    @Override
    public boolean hasTraversalFailed(TraversalState<CortexVertex> s) {
        return false;
    }
}
