package uk.ac.ox.well.indiana.commands.gg;

import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

public abstract class AbstractTraversalStopper implements TraversalStopper {
    public int maxJunctions = 2;
    public int distanceToGoal = Integer.MAX_VALUE;

    public boolean keepGoing(CortexRecord cr, DirectedGraph<String, DefaultEdge> g, int junctions) {
        return !hasTraversalSucceeded(cr, g, junctions) && !hasTraversalFailed(cr, g, junctions);
    }

    public int getDistanceToGoal() { return distanceToGoal; }
}
