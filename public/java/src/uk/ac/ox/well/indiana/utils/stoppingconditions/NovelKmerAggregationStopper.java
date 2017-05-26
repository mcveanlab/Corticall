package uk.ac.ox.well.indiana.utils.stoppingconditions;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedEdge;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedVertex;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;

import java.util.Set;

/**
 * Created by kiran on 24/03/2017.
 */
public class NovelKmerAggregationStopper extends AbstractTraversalStopper<CortexVertex, CortexEdge> {
    private boolean haveSeenNovelKmers = false;

    @Override
    public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<CortexVertex, CortexEdge> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors) {
        return cr == null || cr.getInDegree(childColor) == 0 || cr.getOutDegree(childColor) == 0;
    }

    @Override
    public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<CortexVertex, CortexEdge> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors) {
        return false;
    }

    @Override
    public boolean hasTraversalSucceeded(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedGraph<CortexVertex, CortexEdge> previousGraph, CortexGraph rois) {
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
    public boolean hasTraversalFailed(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedGraph<CortexVertex, CortexEdge> previousGraph, CortexGraph rois) {
        return currentGraphSize >= 100 && !haveSeenNovelKmers;
    }
}
