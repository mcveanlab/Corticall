package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalState;

/**
 * Created by kiran on 20/10/2017.
 */
public class NovelContinuationStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    private boolean startedWithANovelKmer = false;
    private int numKmersSeen = 0;

    @Override
    public boolean hasTraversalSucceeded(TraversalState<CortexVertex> s) {
        if (s.getCurrentJunctionDepth() > 0 && numKmersSeen <= s.getCurrentVertex().getKmerAsString().length() && s.getRois().findRecord(s.getCurrentVertex().getKmerAsByteKmer()) != null) {
            startedWithANovelKmer = true;
        }
        numKmersSeen++;

        return (s.childBranchesAlreadyTraversed() && s.getNumAdjacentEdges() != 1) || s.reachedMaxBranchLength();
    }

    @Override
    public boolean hasTraversalFailed(TraversalState<CortexVertex> s) {
        return (s.getCurrentJunctionDepth() > 0 && !startedWithANovelKmer) || s.getCurrentJunctionDepth() > 3;
    }
}
