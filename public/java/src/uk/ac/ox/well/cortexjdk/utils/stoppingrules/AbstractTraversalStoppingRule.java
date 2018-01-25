package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalState;

public abstract class AbstractTraversalStoppingRule<V, E> implements TraversalStoppingRule<V, E> {
    private boolean traversalSucceeded = false;
    private boolean traversalFailed = false;

    @Override
    public boolean keepGoing(TraversalState<V> s) {
        traversalSucceeded = hasTraversalSucceeded(s);
        traversalFailed = hasTraversalFailed(s);

        return !traversalSucceeded && !traversalFailed;
    }

    @Override
    public boolean hasTraversalSucceeded(TraversalState<V> s) {
        return false;
    }

    @Override
    public boolean hasTraversalFailed(TraversalState<V> s) {
        return true;
    }

    public boolean traversalSucceeded() { return traversalSucceeded; }
    public boolean traversalFailed() { return traversalFailed; }
}
