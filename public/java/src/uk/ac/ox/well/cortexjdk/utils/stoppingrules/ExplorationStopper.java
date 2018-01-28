package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalState;

public class ExplorationStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    private boolean novelKmerFound = false;
    private int distanceFromLastNovelKmer = -1;

    @Override
    public boolean hasTraversalSucceeded(TraversalState<CortexVertex> s) {
        boolean childHasCoverage = s.getCurrentVertex().getCortexRecord().getCoverage(s.getTraversalColor()) > 0;
        boolean parentHasCoverage = false;

        for (int c : s.getJoiningColors()) {
            parentHasCoverage |= s.getCurrentVertex().getCortexRecord().getCoverage(c) > 0;
        }

        if (novelKmerFound) {
            distanceFromLastNovelKmer++;
        }

        if (childHasCoverage && !parentHasCoverage) {
            novelKmerFound = true;
            distanceFromLastNovelKmer = 0;
        }

        return novelKmerFound && (distanceFromLastNovelKmer > 300 || s.childBranchesAlreadyTraversed());
    }

    @Override
    public boolean hasTraversalFailed(TraversalState<CortexVertex> s) {
        return !novelKmerFound && (s.getCurrentGraphSize() > 500 || s.getNumAdjacentEdges() == 0 || s.getCurrentJunctionDepth() > 5);
    }
}
