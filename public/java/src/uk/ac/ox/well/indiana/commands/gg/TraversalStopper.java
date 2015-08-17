package uk.ac.ox.well.indiana.commands.gg;

import org.jgrapht.DirectedGraph;
import org.jgrapht.WeightedGraph;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

public interface TraversalStopper<V, E> {
    boolean keepGoing(CortexRecord cr, DirectedGraph<V, E> g, int junctions, int size);

    boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<V, E> g, int junctions, int size);

    boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<V, E> g, int junctions, int size);

    int maxJunctionsAllowed();
}
