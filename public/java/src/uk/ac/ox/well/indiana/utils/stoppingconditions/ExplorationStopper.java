package uk.ac.ox.well.indiana.utils.stoppingconditions;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedEdge;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedVertex;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;

import java.util.Set;

public class ExplorationStopper extends AbstractTraversalStopper<CortexVertex, CortexEdge> {
    private boolean novelKmerFound = false;
    private int distanceFromLastNovelKmer = -1;

    @Override
    public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<CortexVertex, CortexEdge> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors) {
        return false;
    }

    @Override
    public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<CortexVertex, CortexEdge> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors) {
        return false;
    }

    @Override
    public boolean hasTraversalSucceeded(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedGraph<CortexVertex, CortexEdge> previousGraph, CortexGraph rois) {
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

        return novelKmerFound && (distanceFromLastNovelKmer > 100 || childrenAlreadyTraversed);
    }

    @Override
    public boolean hasTraversalFailed(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedGraph<CortexVertex, CortexEdge> previousGraph, CortexGraph rois) {
        return !novelKmerFound && (currentGraphSize > 200 || numAdjacentEdges == 0 || currentTraversalDepth > 3);
    }
}
