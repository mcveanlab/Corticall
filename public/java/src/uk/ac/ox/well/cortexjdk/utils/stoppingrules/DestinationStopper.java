package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalState;

public class DestinationStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    @Override
    public boolean hasTraversalSucceeded(TraversalState<CortexVertex> s) {
        return s.getSinks().contains(s.getCurrentVertex().getKmerAsString());
    }

    @Override
    public boolean hasTraversalFailed(TraversalState<CortexVertex> s) {
        float coeff = 5.0f;
        int junctionLimit = 1 + (int) Math.ceil(coeff * Math.exp(-0.0001 * (double) s.getCurrentGraphSize()));

        return s.getCurrentJunctionDepth() > junctionLimit || s.reachedMaxBranchLength();
    }
}
