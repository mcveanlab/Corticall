package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;

import java.util.Set;

public class ExplorationStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    private boolean novelKmerFound = false;
    private int distanceFromLastNovelKmer = -1;

    @Override
    public boolean hasTraversalSucceeded(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousGraph, DeBruijnGraph rois) {
        boolean childHasCoverage = cv.getCr().getCoverage(traversalColor) > 0;
        boolean parentHasCoverage = false;

        for (int c : joiningColors) {
            parentHasCoverage |= cv.getCr().getCoverage(c) > 0;
        }

        if (novelKmerFound) {
            distanceFromLastNovelKmer++;
        }

        if (childHasCoverage && !parentHasCoverage) {
            novelKmerFound = true;
            distanceFromLastNovelKmer = 0;
        }

        return novelKmerFound && (distanceFromLastNovelKmer > 300 || childrenAlreadyTraversed);
    }

    @Override
    public boolean hasTraversalFailed(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousGraph, DeBruijnGraph rois) {
        return !novelKmerFound && (currentGraphSize > 500 || numAdjacentEdges == 0 || currentTraversalDepth > 5);
    }
}
