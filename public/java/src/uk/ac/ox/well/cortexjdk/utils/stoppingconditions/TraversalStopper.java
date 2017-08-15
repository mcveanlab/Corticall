package uk.ac.ox.well.cortexjdk.utils.stoppingconditions;

import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;

import java.util.Set;

public interface TraversalStopper<V, E> {
    // This method tells the traversal engine to keep exploring - it hasn't finished evaluating the branch(es) of the graph yet
    boolean keepGoing(V cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<V, E> previousGraph, CortexGraph rois);

    // This method tells the traversal engine to stop and accept the graph branch
    boolean hasTraversalSucceeded(V cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<V, E> previousGraph, CortexGraph rois);

    // This method tells the traversal engine to stop and reject the graph branch
    boolean hasTraversalFailed(V cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<V, E> previousGraph, CortexGraph rois);

    boolean traversalSucceeded();
    boolean traversalFailed();
}
