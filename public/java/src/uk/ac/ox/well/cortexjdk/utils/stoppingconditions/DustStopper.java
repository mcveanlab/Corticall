package uk.ac.ox.well.cortexjdk.utils.stoppingconditions;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.sequence.CortexUtils;
import uk.ac.ox.well.cortexjdk.utils.traversal.AnnotatedEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.AnnotatedVertex;

import java.util.Set;

public class DustStopper extends AbstractTraversalStopper<AnnotatedVertex, AnnotatedEdge> {
    private int sinceLastLowComplexity = 0;

    @Override
    public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors) {
        // We've succeeded if we reach the parents or we run out of edges to traverse
        boolean hasNoIncomingEdges = cr.getInDegree(childColor) == 0;
        boolean hasNoOutgoingEdges = cr.getOutDegree(childColor) == 0;

        boolean reunion = false;
        for (int c : parentColors) {
            reunion |= cr.getCoverage(c) > 0;
        }

        return hasNoIncomingEdges || hasNoOutgoingEdges || reunion;
    }

    @Override
    public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors) {
        // We've failed if we end up seeing a lot of non low-complexity kmers
        if (CortexUtils.isLowComplexity(cr, childColor)) {
            sinceLastLowComplexity = 0;
        } else {
            sinceLastLowComplexity++;
        }

        return sinceLastLowComplexity >= cr.getKmerSize();
    }
}
