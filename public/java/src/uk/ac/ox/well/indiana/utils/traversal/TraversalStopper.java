package uk.ac.ox.well.indiana.utils.traversal;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.util.Set;

public interface TraversalStopper<V, E> {
    boolean keepGoing(CortexRecord cr, DirectedGraph<V, E> g, int junctions, int size, int edges, Set<Integer> childColors, Set<Integer> parentColors);

    boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<V, E> g, int junctions, int size, int edges, Set<Integer> childColors, Set<Integer> parentColors);

    boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<V, E> g, int junctions, int size, int edges, Set<Integer> childColors, Set<Integer> parentColors);

    int maxJunctionsAllowed();

    boolean anyTraversalSucceeded();
    boolean traversalSucceeded();
    boolean traversalFailed();
}
