package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalState;

public class ContaminantStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    @Override
    public boolean hasTraversalSucceeded(TraversalState<CortexVertex> s) {
        boolean parentsHaveCoverage = false;

        for (int c : s.getJoiningColors()) {
            parentsHaveCoverage |= s.getCurrentVertex().getCortexRecord().getCoverage(c) > 0;
        }

        return s.getCurrentVertex().getCanonicalKmer() != null && (parentsHaveCoverage || s.getNumAdjacentEdges() == 0);
    }

    @Override
    public boolean hasTraversalFailed(TraversalState<CortexVertex> s) {
        boolean parentsHaveCoverage = false;

        for (int c : s.getJoiningColors()) {
            parentsHaveCoverage |= s.getCurrentVertex().getCortexRecord().getCoverage(c) > 0;
        }

        return s.getCurrentVertex().getCortexRecord() != null && parentsHaveCoverage;
    }
}
