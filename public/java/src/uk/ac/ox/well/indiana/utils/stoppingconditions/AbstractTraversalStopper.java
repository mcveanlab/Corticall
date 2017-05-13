package uk.ac.ox.well.indiana.utils.stoppingconditions;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;

import java.util.Set;

public abstract class AbstractTraversalStopper<V, E> implements TraversalStopper<V, E> {
    private boolean anyTraversalSucceeded = false;
    private boolean traversalSucceeded = false;
    private boolean traversalFailed = false;

    public boolean keepGoing(CortexRecord cr, DirectedGraph<V, E> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors) {
        traversalSucceeded = hasTraversalSucceeded(cr, g, junctions, size, edges, childColor, parentColors);
        traversalFailed = hasTraversalFailed(cr, g, junctions, size, edges, childColor, parentColors);

        if (traversalSucceeded) {
            anyTraversalSucceeded = true;
        }

        return !traversalSucceeded && !traversalFailed;
    }

    @Override
    public boolean keepGoing(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, DirectedGraph<V, E> previousGraph) {
        traversalSucceeded = hasTraversalSucceeded(cv, goForward, traversalColor, joiningColors, currentTraversalDepth, currentGraphSize, numAdjacentEdges, previousGraph);
        traversalFailed = hasTraversalFailed(cv, goForward, traversalColor, joiningColors, currentTraversalDepth, currentGraphSize, numAdjacentEdges, previousGraph);

        if (traversalSucceeded) {
            anyTraversalSucceeded = true;
        }

        return !traversalSucceeded && !traversalFailed;
    }

    @Override
    public boolean hasTraversalSucceeded(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, DirectedGraph<V, E> previousGraph) {
        return false;
    }

    @Override
    public boolean hasTraversalFailed(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, DirectedGraph<V, E> previousGraph) {
        return false;
    }

    public boolean anyTraversalSucceeded() { return anyTraversalSucceeded; }
    public boolean traversalSucceeded() { return traversalSucceeded; }
    public boolean traversalFailed() { return traversalFailed; }
}
