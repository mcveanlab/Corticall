package uk.ac.ox.well.indiana.utils.stoppingconditions;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedEdge;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedVertex;

import java.util.Set;

/**
 * Created by kiran on 08/05/2017.
 */
public class NahrStopper extends AbstractTraversalStopper<AnnotatedVertex, AnnotatedEdge> {
    int distanceFromGoal = -1;

    @Override
    public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors) {
        boolean childHasCoverage = cr.getCoverage(childColor) > 0;
        boolean aParentHasCoverage = false;

        for (int c : parentColors) {
            aParentHasCoverage |= cr.getCoverage(c) > 0;
        }

        if (childHasCoverage && !aParentHasCoverage) {
            distanceFromGoal = 0;
        }

        return distanceFromGoal > 0 && (cr.getInDegree(childColor) == 0 || cr.getOutDegree(childColor) == 0 || distanceFromGoal > 500);

        //return false;
    }

    @Override
    public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors) {
        return false;
    }
}
