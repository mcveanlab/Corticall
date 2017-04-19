package uk.ac.ox.well.indiana.utils.traversal;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.util.Set;

/**
 * Created by kiran on 20/01/2016.
 */
public class ContaminantStopper extends AbstractTraversalStopper<AnnotatedVertex, AnnotatedEdge> {
    public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int depth, int size, int edges, int childColor, Set<Integer> parentColors) {
        return cr != null && (cr.getCoverage(1) > 0 || cr.getCoverage(2) > 0 || edges == 0);
    }

    @Override
    public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int depth, int size, int edges, int childColor, Set<Integer> parentColors) {
        return cr != null && (cr.getCoverage(1) > 0 || cr.getCoverage(2) > 0);
    }

    @Override
    public int maxJunctionsAllowed() {
        return 0;
    }
}
