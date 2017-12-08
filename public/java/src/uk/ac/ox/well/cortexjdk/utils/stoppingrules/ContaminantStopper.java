package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;

import java.util.Set;

public class ContaminantStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    @Override
    public boolean hasTraversalSucceeded(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousGraph, DeBruijnGraph rois) {
        boolean parentsHaveCoverage = false;

        for (int c : joiningColors) {
            parentsHaveCoverage |= cv.getCortexRecord().getCoverage(c) > 0;
        }

        return cv.getCanonicalKmer() != null && (parentsHaveCoverage || numAdjacentEdges == 0);
    }

    @Override
    public boolean hasTraversalFailed(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousGraph, DeBruijnGraph rois) {
        boolean parentsHaveCoverage = false;

        for (int c : joiningColors) {
            parentsHaveCoverage |= cv.getCortexRecord().getCoverage(c) > 0;
        }

        return cv.getCortexRecord() != null && parentsHaveCoverage;
    }
}
