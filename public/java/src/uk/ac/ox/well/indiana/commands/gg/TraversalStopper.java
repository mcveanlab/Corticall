package uk.ac.ox.well.indiana.commands.gg;

import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

/**
 * Created by kiran on 21/06/2015.
 */
public interface TraversalStopper {
    boolean keepGoing(CortexRecord cr, DirectedGraph<String, DefaultEdge> g, int depth);
}
