package uk.ac.ox.well.cortexjdk.utils.stoppingconditions;

import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;

import java.util.Set;

public class DustStopper extends AbstractTraversalStopper<CortexVertex, CortexEdge> {
    private int sinceLastLowComplexity = 0;

    @Override
    public boolean hasTraversalSucceeded(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousGraph, DeBruijnGraph rois) {
        // We've succeeded if we reach the parents or we run out of edges to traverse
        boolean hasNoIncomingEdges = cv.getCr().getInDegree(traversalColor) == 0;
        boolean hasNoOutgoingEdges = cv.getCr().getOutDegree(traversalColor) == 0;

        boolean reunion = false;
        for (int c : joiningColors) {
            reunion |= cv.getCr().getCoverage(c) > 0;
        }

        return hasNoIncomingEdges || hasNoOutgoingEdges || reunion;
    }

    @Override
    public boolean hasTraversalFailed(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousGraph, DeBruijnGraph rois) {
        // We've failed if we end up seeing a lot of non low-complexity kmers
        if (isLowComplexity(cv.getCr(), traversalColor)) {
            sinceLastLowComplexity = 0;
        } else {
            sinceLastLowComplexity++;
        }

        return sinceLastLowComplexity >= cv.getCr().getKmerSize();
    }

    private boolean isLowComplexity(CortexRecord cr, int traversalColor) {
        return cr.getInDegree(traversalColor) + cr.getOutDegree(traversalColor) > 4;
    }
}
