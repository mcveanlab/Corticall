package uk.ac.ox.well.indiana.utils.traversal;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

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

    public boolean anyTraversalSucceeded() { return anyTraversalSucceeded; }
    public boolean traversalSucceeded() { return traversalSucceeded; }
    public boolean traversalFailed() { return traversalFailed; }
}
