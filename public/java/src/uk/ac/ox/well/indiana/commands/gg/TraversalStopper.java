package uk.ac.ox.well.indiana.commands.gg;

import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

public interface TraversalStopper<V, E> {
    boolean keepGoing(CortexRecord cr, DirectedGraph<V, E> g, int junctions);

    boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<V, E> g, int junctions);

    boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<V, E> g, int junctions);

    int maxJunctionsAllowed();

    int getDistanceToGoal();
}
