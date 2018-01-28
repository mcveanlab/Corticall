package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalState;

/**
 * Created by kiran on 08/05/2017.
 */
public class NahrStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    private boolean foundNovels = false;
    private int distanceFromLastNovel = 0;

    @Override
    public boolean hasTraversalSucceeded(TraversalState<CortexVertex> s) {
        if (foundNovels) {
            distanceFromLastNovel++;
        }

        if (s.getRois().findRecord(s.getCurrentVertex().getCortexRecord().getCanonicalKmer()) != null) {
            foundNovels = true;
            distanceFromLastNovel++;
        }

        return foundNovels && (distanceFromLastNovel >= 1000 || s.getCurrentJunctionDepth() >= 5 || s.getNumAdjacentEdges() == 0 || s.childBranchesAlreadyTraversed());
    }

    @Override
    public boolean hasTraversalFailed(TraversalState<CortexVertex> s) {
        return !foundNovels && (s.getCurrentGraphSize() >= 1000 || s.getCurrentJunctionDepth() >= 2 || s.getNumAdjacentEdges() == 0);
    }
}
