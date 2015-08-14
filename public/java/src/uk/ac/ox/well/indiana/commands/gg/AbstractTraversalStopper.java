package uk.ac.ox.well.indiana.commands.gg;

import org.jgrapht.DirectedGraph;
import org.jgrapht.WeightedGraph;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

public abstract class AbstractTraversalStopper<V, E> implements TraversalStopper<V, E> {
    public boolean keepGoing(CortexRecord cr, DefaultDirectedWeightedGraph<V, E> g, int junctions, int size) {
        return !hasTraversalSucceeded(cr, g, junctions, size) && !hasTraversalFailed(cr, g, junctions, size);
    }
}
