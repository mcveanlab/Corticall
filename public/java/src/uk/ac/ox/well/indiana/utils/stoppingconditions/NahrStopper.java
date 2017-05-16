package uk.ac.ox.well.indiana.utils.stoppingconditions;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;

import java.util.Set;

/**
 * Created by kiran on 08/05/2017.
 */
public class NahrStopper extends AbstractTraversalStopper<CortexVertex, CortexEdge> {
    private boolean seenNovelKmer = false;
    private int distanceSinceLastNovel = 0;

    @Override
    public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<CortexVertex, CortexEdge> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors) {
        return false;
    }

    @Override
    public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<CortexVertex, CortexEdge> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors) {
        return false;
    }

    @Override
    public boolean hasTraversalSucceeded(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, DirectedGraph<CortexVertex, CortexEdge> previousGraph, CortexGraph rois) {
        if (seenNovelKmer) {
            distanceSinceLastNovel++;
        }

        if (rois.findRecord(cv.getCr().getCortexKmer()) != null) {
            seenNovelKmer = true;
            distanceSinceLastNovel = 0;
        }

        return seenNovelKmer && (distanceSinceLastNovel > 1000 || numAdjacentEdges == 0);
    }

    @Override
    public boolean hasTraversalFailed(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, DirectedGraph<CortexVertex, CortexEdge> previousGraph, CortexGraph rois) {
        return !seenNovelKmer && (currentTraversalDepth > 1000 || numAdjacentEdges == 0);
    }
}
