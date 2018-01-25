package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalState;

import java.util.Set;

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
