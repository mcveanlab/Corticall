package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;

import java.util.Set;

/**
 * Created by kiran on 20/10/2017.
 */
public class NovelContinuationStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    private int novelKmersSeen = 0;
    private int distanceSinceJoin = 0;
    private boolean hasJoined = false;

    @Override
    public boolean hasTraversalSucceeded(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousGraph, DeBruijnGraph rois) {
        if (rois.findRecord(cv.getBk()) != null) {
            novelKmersSeen++;
        }

        /*
        if (hasJoined) {
            distanceSinceJoin++;
        }

        for (int c : joiningColors) {
            hasJoined |= (cv.getCr().getCoverage(c) > 0);
        }
        */

        return novelKmersSeen > 0 && (childrenAlreadyTraversed ? numAdjacentEdges > 1 : numAdjacentEdges == 0);

        //return novelKmersSeen > 0 && hasJoined && (distanceSinceJoin >= 30 || numAdjacentEdges != 1);
    }

    @Override
    public boolean hasTraversalFailed(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousGraph, DeBruijnGraph rois) {
        return novelKmersSeen == 0 && (currentTraversalDepth >= 3 || numAdjacentEdges == 0);
        //return novelKmersSeen == 0 && (currentTraversalDepth > 1 || numAdjacentEdges == 0);
    }
}
