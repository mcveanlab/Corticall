package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;

import java.util.Set;

/**
 * Created by kiran on 24/03/2017.
 */
public class NovelKmerAggregationStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    private boolean haveSeenNovelKmers = false;

    @Override
    public boolean hasTraversalSucceeded(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousGraph, DeBruijnGraph rois) {
        boolean childHasCoverage = cv.getCr().getCoverage(traversalColor) > 0;
        boolean parentsHaveCoverage = false;

        for (int c : joiningColors) {
            parentsHaveCoverage |= cv.getCr().getCoverage(c) > 0;
        }

        if (childHasCoverage && !parentsHaveCoverage) {
            haveSeenNovelKmers = true;
        }

        return haveSeenNovelKmers && parentsHaveCoverage;
    }

    @Override
    public boolean hasTraversalFailed(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousGraph, DeBruijnGraph rois) {
        return !haveSeenNovelKmers && (currentGraphSize >= 100 || currentTraversalDepth >= 3);
    }
}
