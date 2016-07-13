package uk.ac.ox.well.indiana.utils.traversal;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

public interface TraversalStopper<V, E> {
    boolean keepGoing(CortexRecord cr, DirectedGraph<V, E> g, int junctions, int size, int edges);

    boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<V, E> g, int junctions, int size, int edges);

    boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<V, E> g, int junctions, int size, int edges);

    int maxJunctionsAllowed();

    boolean anyTraversalSucceeded();
    boolean traversalSucceeded();
    boolean traversalFailed();
}
