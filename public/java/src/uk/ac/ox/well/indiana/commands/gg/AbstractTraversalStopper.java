package uk.ac.ox.well.indiana.commands.gg;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

public abstract class AbstractTraversalStopper<V, E> implements TraversalStopper<V, E> {
    private boolean anyTraversalSucceeded = false;
    private boolean traversalSucceeded = false;
    private boolean traversalFailed = false;

    public boolean keepGoing(CortexRecord cr, DirectedGraph<V, E> g, int junctions, int size, int edges) {
        //return !hasTraversalSucceeded(cr, g, junctions, size) && !hasTraversalFailed(cr, g, junctions, size);
        traversalSucceeded = hasTraversalSucceeded(cr, g, junctions, size, edges);
        traversalFailed = hasTraversalFailed(cr, g, junctions, size, edges);

        if (traversalSucceeded) {
            anyTraversalSucceeded = true;
        }

        return !traversalSucceeded && !traversalFailed;
    }

    public boolean anyTraversalSucceeded() { return anyTraversalSucceeded; }
    public boolean traversalSucceeded() { return traversalSucceeded; }
    public boolean traversalFailed() { return traversalFailed; }
}
