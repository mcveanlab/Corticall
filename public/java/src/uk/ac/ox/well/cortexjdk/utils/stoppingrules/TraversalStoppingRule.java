package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalState;

public interface TraversalStoppingRule<V, E> {
    // This method tells the traversal engine to keep exploring - it hasn't finished evaluating the branch(es) of the graph yet
    boolean keepGoing(TraversalState<V> s);

    // This method tells the traversal engine to stop and accept the graph branch
    boolean hasTraversalSucceeded(TraversalState<V> s);

    // This method tells the traversal engine to stop and reject the graph branch
    boolean hasTraversalFailed(TraversalState<V> s);

    boolean traversalSucceeded();
    boolean traversalFailed();
}
