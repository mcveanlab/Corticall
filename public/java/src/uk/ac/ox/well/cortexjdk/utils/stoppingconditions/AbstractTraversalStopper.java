package uk.ac.ox.well.cortexjdk.utils.stoppingconditions;

import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.DeBruijnGraph;

import java.util.Set;

public abstract class AbstractTraversalStopper<V, E> implements TraversalStopper<V, E> {
    private boolean traversalSucceeded = false;
    private boolean traversalFailed = false;

    @Override
    public boolean keepGoing(V cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<V, E> previousGraph, DeBruijnGraph rois) {
        traversalSucceeded = hasTraversalSucceeded(cv, goForward, traversalColor, joiningColors, currentTraversalDepth, currentGraphSize, numAdjacentEdges, childrenAlreadyTraversed, previousGraph, rois);
        traversalFailed = hasTraversalFailed(cv, goForward, traversalColor, joiningColors, currentTraversalDepth, currentGraphSize, numAdjacentEdges, childrenAlreadyTraversed, previousGraph, rois);

        return !traversalSucceeded && !traversalFailed;
    }

    @Override
    public boolean hasTraversalSucceeded(V cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<V, E> previousGraph, DeBruijnGraph rois) {
        return false;
    }

    @Override
    public boolean hasTraversalFailed(V cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<V, E> previousGraph, DeBruijnGraph rois) {
        return true;
    }

    public boolean traversalSucceeded() { return traversalSucceeded; }
    public boolean traversalFailed() { return traversalFailed; }
}
