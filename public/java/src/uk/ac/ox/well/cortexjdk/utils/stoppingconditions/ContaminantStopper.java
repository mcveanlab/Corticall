package uk.ac.ox.well.cortexjdk.utils.stoppingconditions;

import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;

import java.util.Set;

public class ContaminantStopper extends AbstractTraversalStopper<CortexVertex, CortexEdge> {
    @Override
    public boolean hasTraversalSucceeded(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousGraph, CortexGraph rois) {
        boolean parentsHaveCoverage = false;

        for (int c : joiningColors) {
            parentsHaveCoverage |= cv.getCr().getCoverage(c) > 0;
        }

        return cv.getCk() != null && (parentsHaveCoverage || numAdjacentEdges == 0);
    }

    @Override
    public boolean hasTraversalFailed(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousGraph, CortexGraph rois) {
        boolean parentsHaveCoverage = false;

        for (int c : joiningColors) {
            parentsHaveCoverage |= cv.getCr().getCoverage(c) > 0;
        }

        return cv.getCr() != null && parentsHaveCoverage;
    }
}
