package uk.ac.ox.well.indiana.commands.gg;

import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

public abstract class AbstractTraversalStopper implements TraversalStopper {
    public boolean keepGoing(CortexRecord cr, DirectedGraph<String, DefaultEdge> g, int depth) {
        return !hasTraversalSucceeded(cr, g, depth) && !hasTraversalFailed(cr, g, depth);
    }
}
