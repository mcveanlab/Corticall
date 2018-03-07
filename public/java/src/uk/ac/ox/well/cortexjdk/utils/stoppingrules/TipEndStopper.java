package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalState;

public class TipEndStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    @Override
    public boolean hasTraversalSucceeded(TraversalState<CortexVertex> s) {
        // We should accept this branch if we run out of edges to navigate

        boolean hasNoIncomingEdges = false;
        boolean hasNoOutgoingEdges = false;

        for (int c : s.getTraversalColors()) {
            hasNoIncomingEdges |= s.getCurrentVertex().getCortexRecord().getInDegree(c) == 0;
            hasNoOutgoingEdges |= s.getCurrentVertex().getCortexRecord().getOutDegree(c) == 0;
        }

        return hasNoIncomingEdges || hasNoOutgoingEdges;
    }

    @Override
    public boolean hasTraversalFailed(TraversalState<CortexVertex> s) {
        // We should reject this branch if we make it to the end of our traversal and we reconnect with a parent

        boolean reunion = false;
        for (int c : s.getJoiningColors()) {
            reunion |= s.getCurrentVertex().getCortexRecord().getCoverage(c) > 0;
        }

        return reunion;
    }
}
