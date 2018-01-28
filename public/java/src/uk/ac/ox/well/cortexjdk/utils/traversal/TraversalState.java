package uk.ac.ox.well.cortexjdk.utils.traversal;

import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.DeBruijnGraph;

import java.util.Set;

final public class TraversalState<V> {
    final private V currentVertex;
    final private boolean goForward;
    final private int traversalColor;
    final private Set<Integer> joiningColors;
    final private int currentJunctionDepth;
    final private int currentGraphSize;
    final private int numAdjacentEdges;
    final private boolean childrenAlreadyTraversed;
    final private DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousGraph;
    final private DeBruijnGraph rois;
    final private boolean reachedMaxBranchLength;

    public TraversalState(V currentVertex,
                          boolean goForward,
                          int traversalColor,
                          Set<Integer> joiningColors,
                          int currentJunctionDepth,
                          int currentGraphSize,
                          int numAdjacentEdges,
                          boolean childrenAlreadyTraversed,
                          boolean reachedMaxBranchLength,
                          DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousGraph,
                          DeBruijnGraph rois) {
        this.currentVertex = currentVertex;
        this.goForward = goForward;
        this.traversalColor = traversalColor;
        this.joiningColors = joiningColors;
        this.currentJunctionDepth = currentJunctionDepth;
        this.currentGraphSize = currentGraphSize;
        this.numAdjacentEdges = numAdjacentEdges;
        this.childrenAlreadyTraversed = childrenAlreadyTraversed;
        this.previousGraph = previousGraph;
        this.rois = rois;
        this.reachedMaxBranchLength = reachedMaxBranchLength;
    }

    public V getCurrentVertex() { return currentVertex; }

    public boolean goForward() { return goForward; }

    public int getTraversalColor() { return traversalColor; }

    public Set<Integer> getJoiningColors() { return joiningColors; }

    public int getCurrentJunctionDepth() { return currentJunctionDepth; }

    public int getCurrentGraphSize() { return currentGraphSize; }

    public int getNumAdjacentEdges() { return numAdjacentEdges; }

    public boolean childBranchesAlreadyTraversed() { return childrenAlreadyTraversed; }

    public boolean reachedMaxBranchLength() { return reachedMaxBranchLength; }

    public DirectedWeightedPseudograph<CortexVertex, CortexEdge> getPreviousGraph() { return previousGraph; }

    public DeBruijnGraph getRois() { return rois; }
}
