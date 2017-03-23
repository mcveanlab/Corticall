package uk.ac.ox.well.indiana.utils.traversal;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.util.Set;

public class ExplorationStopper extends AbstractTraversalStopper<AnnotatedVertex, AnnotatedEdge> {
    @Override
    public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions, int size, int edges, Set<Integer> childColors, Set<Integer> parentColors) {
        return size >= 50 || junctions >= 1 || edges == 0;
    }

    @Override
    public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions, int size, int edges, Set<Integer> childColors, Set<Integer> parentColors) {
        return false;
    }

    @Override
    public int maxJunctionsAllowed() {
        return 0;
    }
}
