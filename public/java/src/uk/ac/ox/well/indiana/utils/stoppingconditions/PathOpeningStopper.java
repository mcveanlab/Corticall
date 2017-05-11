package uk.ac.ox.well.indiana.utils.stoppingconditions;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedEdge;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedVertex;

import java.util.Set;

/**
 * Created by kiran on 07/05/2017.
 */
public class PathOpeningStopper extends AbstractTraversalStopper<AnnotatedVertex, AnnotatedEdge>{
    @Override
    public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors) {
        return (cr.getInDegree(childColor) > 1 || cr.getOutDegree(childColor) > 1);
    }

    @Override
    public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors) {
        return (cr.getInDegree(childColor) == 0 || cr.getOutDegree(childColor) == 0);
    }
}
