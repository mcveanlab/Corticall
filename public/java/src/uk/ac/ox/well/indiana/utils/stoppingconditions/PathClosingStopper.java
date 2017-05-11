package uk.ac.ox.well.indiana.utils.stoppingconditions;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedEdge;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedVertex;

import java.util.Set;

/**
 * Created by kiran on 08/05/2017.
 */
public class PathClosingStopper extends AbstractTraversalStopper<AnnotatedVertex, AnnotatedEdge> {
    @Override
    public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors) {
        String fw = cr.getKmerAsString();
        String rc = SequenceUtils.reverseComplement(fw);

        boolean rejoined = g.containsVertex(new AnnotatedVertex(fw)) || g.containsVertex(new AnnotatedVertex(rc));

        return size >= 2 && rejoined;
    }

    @Override
    public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors) {
        return size >= 100;
    }
}
