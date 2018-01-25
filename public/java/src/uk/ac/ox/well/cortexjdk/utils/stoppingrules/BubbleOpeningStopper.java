package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalState;

import java.util.Set;

/**
 * Created by kiran on 07/05/2017.
 */
public class BubbleOpeningStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    private int novelKmersSeen = 0;
    private int distanceSinceJoin = 0;
    private boolean hasJoined = false;

    @Override
    public boolean hasTraversalSucceeded(TraversalState<CortexVertex> s) {
        if (s.getRois().findRecord(s.getCurrentVertex().getKmerAsByteKmer()) != null) {
            novelKmersSeen++;
        }

        if (hasJoined) {
            distanceSinceJoin++;
        }

        for (int c : s.getJoiningColors()) {
            hasJoined |= (s.getCurrentVertex().getCortexRecord().getCoverage(c) > 0);
        }

        return novelKmersSeen > 0 && hasJoined && (distanceSinceJoin >= 30 || s.getNumAdjacentEdges() != 1);
    }

    @Override
    public boolean hasTraversalFailed(TraversalState<CortexVertex> s) {
        return novelKmersSeen == 0 && (s.getCurrentTraversalDepth() >= 5 || s.getNumAdjacentEdges() == 0);
    }
}
