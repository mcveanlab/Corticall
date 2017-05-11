package uk.ac.ox.well.indiana.utils.stoppingconditions;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;

import java.util.Set;

public interface TraversalStopper<V, E> {
    // This method tells the traversal engine to keep exploring - it hasn't finished evaluating the branch(es) of the graph yet
    boolean keepGoing(CortexRecord cr, DirectedGraph<V, E> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors);
    boolean keepGoing(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, DirectedGraph<V, E> previousGraph);

    // This method tells the traversal engine to stop and accept the graph branch
    boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<V, E> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors);
    boolean hasTraversalSucceeded(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, DirectedGraph<V, E> previousGraph);

    // This method tells the traversal engine to stop and reject the graph branch
    boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<V, E> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors);
    boolean hasTraversalFailed(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, DirectedGraph<V, E> previousGraph);

    boolean anyTraversalSucceeded();
    boolean traversalSucceeded();
    boolean traversalFailed();
}
