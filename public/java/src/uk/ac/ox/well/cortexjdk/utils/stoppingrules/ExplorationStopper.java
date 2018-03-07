package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalState;

public class ExplorationStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    @Override
    public boolean hasTraversalSucceeded(TraversalState<CortexVertex> s) {
        return (s.reachedMaxBranchLength() || s.getNumAdjacentEdges() == 0 || s.getCurrentJunctionDepth() >= 3);
    }

    @Override
    public boolean hasTraversalFailed(TraversalState<CortexVertex> s) {
        return false;
    }
}
