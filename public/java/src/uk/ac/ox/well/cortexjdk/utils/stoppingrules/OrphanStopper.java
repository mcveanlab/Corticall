package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;

import java.util.Set;

public class OrphanStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    @Override
    public boolean hasTraversalSucceeded(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousGraph, DeBruijnGraph rois) {
        // We should accept this branch if we make it to the end of our traversal and there are no more edges to navigate

        boolean hasNoIncomingEdges = cv.getCr().getInDegree(traversalColor) == 0;
        boolean hasNoOutgoingEdges = cv.getCr().getOutDegree(traversalColor) == 0;

        return hasNoIncomingEdges || hasNoOutgoingEdges;
    }

    @Override
    public boolean hasTraversalFailed(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousGraph, DeBruijnGraph rois) {
        // We should reject this branch if we ever reconnect with the parental colors.

        boolean reunion = false;
        for (int c : joiningColors) {
            reunion |= cv.getCr().getCoverage(c) > 0;
        }

        return reunion;
    }
}
