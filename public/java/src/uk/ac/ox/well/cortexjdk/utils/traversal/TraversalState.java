package uk.ac.ox.well.cortexjdk.utils.traversal;

import uk.ac.ox.well.cortexjdk.utils.io.graph.DeBruijnGraph;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

final public class TraversalState<V> {
    final private V currentVertex;
    final private boolean goForward;
    final private Set<Integer> traversalColors;
    final private Set<Integer> joiningColors;
    final private int currentJunctionDepth;
    final private int currentGraphSize;
    final private int currentBranchSize;
    final private int numAdjacentEdges;
    final private boolean childrenAlreadyTraversed;
    final private DeBruijnGraph rois;
    final private boolean reachedMaxBranchLength;
    final private Set<String> sinks;

    public TraversalState(V currentVertex,
                          boolean goForward,
                          Set<Integer> traversalColors,
                          Set<Integer> joiningColors,
                          int currentGraphSize,
                          int currentJunctionDepth,
                          int currentBranchSize,
                          int numAdjacentEdges,
                          boolean childrenAlreadyTraversed,
                          boolean reachedMaxBranchLength,
                          DeBruijnGraph rois,
                          String... sinks
    ) {
        this.currentVertex = currentVertex;
        this.goForward = goForward;
        this.traversalColors = traversalColors;
        this.joiningColors = joiningColors;
        this.currentGraphSize = currentGraphSize;
        this.currentJunctionDepth = currentJunctionDepth;
        this.currentBranchSize = currentBranchSize;
        this.numAdjacentEdges = numAdjacentEdges;
        this.childrenAlreadyTraversed = childrenAlreadyTraversed;
        this.rois = rois;
        this.reachedMaxBranchLength = reachedMaxBranchLength;

        Set<String> s = new HashSet<>();
        s.addAll(Arrays.asList(sinks));
        this.sinks = s;
    }

    public V getCurrentVertex() { return currentVertex; }

    public boolean goForward() { return goForward; }

    public Set<Integer> getTraversalColors() { return traversalColors; }

    public Set<Integer> getJoiningColors() { return joiningColors; }

    public int getCurrentGraphSize() { return currentGraphSize; }

    public int getCurrentJunctionDepth() { return currentJunctionDepth; }

    public int getCurrentBranchSize() { return currentBranchSize; }

    public int getNumAdjacentEdges() { return numAdjacentEdges; }

    public boolean childBranchesAlreadyTraversed() { return childrenAlreadyTraversed; }

    public boolean reachedMaxBranchLength() { return reachedMaxBranchLength; }

    public Set<String> getSinks() { return sinks; }

    public DeBruijnGraph getRois() { return rois; }
}
