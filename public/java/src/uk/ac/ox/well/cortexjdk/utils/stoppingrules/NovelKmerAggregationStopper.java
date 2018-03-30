package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalState;

/**
 * Created by kiran on 24/03/2017.
 */
public class NovelKmerAggregationStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    private boolean haveSeenNovelKmers = false;

    @Override
    public boolean hasTraversalSucceeded(TraversalState<CortexVertex> s) {
        boolean childHasCoverage = false;

        for (int c : s.getTraversalColors()) {
            childHasCoverage |= s.getCurrentVertex().getCortexRecord().getCoverage(c) > 0;
        }

        boolean parentsHaveCoverage = false;

        for (int c : s.getJoiningColors()) {
            parentsHaveCoverage |= s.getCurrentVertex().getCortexRecord().getCoverage(c) > 0;
        }

        if (childHasCoverage && !parentsHaveCoverage) {
            haveSeenNovelKmers = true;
        }

        return haveSeenNovelKmers && parentsHaveCoverage;
    }

    @Override
    public boolean hasTraversalFailed(TraversalState<CortexVertex> s) {
        return !haveSeenNovelKmers && (s.getCurrentBranchSize() >= 100 || s.getCurrentJunctionDepth() >= 3);
    }
}
