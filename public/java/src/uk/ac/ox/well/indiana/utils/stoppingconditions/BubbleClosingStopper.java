package uk.ac.ox.well.indiana.utils.stoppingconditions;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedEdge;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedVertex;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;

import java.util.Set;

/**
 * Created by kiran on 08/05/2017.
 */
public class BubbleClosingStopper extends AbstractTraversalStopper<CortexVertex, CortexEdge> {
    @Override
    public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<CortexVertex, CortexEdge> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors) {
        return false;
    }

    @Override
    public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<CortexVertex, CortexEdge> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors) {
        return false;
    }

    @Override
    public boolean hasTraversalSucceeded(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedGraph<CortexVertex, CortexEdge> previousGraph, CortexGraph rois) {
        return previousGraph.containsVertex(cv);
    }

    @Override
    public boolean hasTraversalFailed(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedGraph<CortexVertex, CortexEdge> previousGraph, CortexGraph rois) {
        return currentGraphSize > 200 || currentTraversalDepth >= 5 || numAdjacentEdges == 0;
    }
}
