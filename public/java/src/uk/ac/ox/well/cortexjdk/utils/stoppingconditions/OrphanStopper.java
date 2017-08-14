package uk.ac.ox.well.cortexjdk.utils.stoppingconditions;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.traversal.AnnotatedEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.AnnotatedVertex;

import java.util.Set;

public class OrphanStopper extends AbstractTraversalStopper<AnnotatedVertex, AnnotatedEdge> {
    public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int depth, int size, int edges, int childColor, Set<Integer> parentColors) {
        // We should accept this branch if we make it to the end of our traversal and there are no more edges to navigate

        boolean hasNoIncomingEdges = cr.getInDegree(childColor) == 0;
        boolean hasNoOutgoingEdges = cr.getOutDegree(childColor) == 0;

        return hasNoIncomingEdges || hasNoOutgoingEdges;
    }

    @Override
    public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int depth, int size, int edges, int childColor, Set<Integer> parentColors) {
        // We should reject this branch if we ever reconnect with the parental colors.

        boolean reunion = false;
        for (int c : parentColors) {
            reunion |= cr.getCoverage(c) > 0;
        }

        return reunion;
    }
}
