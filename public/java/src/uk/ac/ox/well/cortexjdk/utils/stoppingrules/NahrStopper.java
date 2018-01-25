package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalState;

import java.util.Set;

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

        return foundNovels && (distanceFromLastNovel >= 1000 || s.getCurrentTraversalDepth() >= 5 || s.getNumAdjacentEdges() == 0 || s.isChildrenAlreadyTraversed());
    }

    @Override
    public boolean hasTraversalFailed(TraversalState<CortexVertex> s) {
        return !foundNovels && (s.getCurrentGraphSize() >= 1000 || s.getCurrentTraversalDepth() >= 2 || s.getNumAdjacentEdges() == 0);
    }
}
