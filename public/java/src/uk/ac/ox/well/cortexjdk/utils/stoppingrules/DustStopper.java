package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalState;

public class DustStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    private int sinceLastLowComplexity = 0;

    @Override
    public boolean hasTraversalSucceeded(TraversalState<CortexVertex> s) {
        // We've succeeded if we reach the parents or we run out of edges to traverse
        boolean hasNoIncomingEdges = s.getCurrentVertex().getCortexRecord().getInDegree(s.getTraversalColor()) == 0;
        boolean hasNoOutgoingEdges = s.getCurrentVertex().getCortexRecord().getOutDegree(s.getTraversalColor()) == 0;

        boolean reunion = false;
        for (int c : s.getJoiningColors()) {
            reunion |= s.getCurrentVertex().getCortexRecord().getCoverage(c) > 0;
        }

        return hasNoIncomingEdges || hasNoOutgoingEdges || reunion;
    }

    @Override
    public boolean hasTraversalFailed(TraversalState<CortexVertex> s) {
        // We've failed if we end up seeing a lot of non low-complexity kmers
        if (isLowComplexity(s.getCurrentVertex().getCortexRecord(), s.getTraversalColor())) {
            sinceLastLowComplexity = 0;
        } else {
            sinceLastLowComplexity++;
        }

        return sinceLastLowComplexity >= s.getCurrentVertex().getCortexRecord().getKmerSize();
    }

    private boolean isLowComplexity(CortexRecord cr, int traversalColor) {
        return cr.getInDegree(traversalColor) + cr.getOutDegree(traversalColor) > 4;
    }
}
