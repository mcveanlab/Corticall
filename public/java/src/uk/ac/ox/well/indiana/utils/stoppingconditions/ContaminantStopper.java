package uk.ac.ox.well.indiana.utils.stoppingconditions;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedEdge;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedVertex;

import java.util.Set;

public class ContaminantStopper extends AbstractTraversalStopper<AnnotatedVertex, AnnotatedEdge> {
    public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int depth, int size, int edges, int childColor, Set<Integer> parentColors) {
        boolean parentsHaveCoverage = false;

        for (int c : parentColors) {
            parentsHaveCoverage |= cr.getCoverage(c) > 0;
        }

        return cr != null && (parentsHaveCoverage || edges == 0);
    }

    @Override
    public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int depth, int size, int edges, int childColor, Set<Integer> parentColors) {
        boolean parentsHaveCoverage = false;

        for (int c : parentColors) {
            parentsHaveCoverage |= cr.getCoverage(c) > 0;
        }

        return cr != null && parentsHaveCoverage;
    }
}
