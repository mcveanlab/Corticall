package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalState;

/**
 * Created by kiran on 10/05/2017.
 */
public class ContigStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    @Override
    public boolean hasTraversalSucceeded(TraversalState<CortexVertex> s) {
        return s.getNumAdjacentEdges() != 1;
    }

    @Override
    public boolean hasTraversalFailed(TraversalState<CortexVertex> s) {
        return false;
    }
}
