package uk.ac.ox.well.cortexjdk.utils.stoppingconditions;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * Created by kiran on 08/05/2017.
 */
public class ShoreStopper extends AbstractTraversalStopper<CortexVertex, CortexEdge> {
    @Override
    public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<CortexVertex, CortexEdge> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors) {
        return false;
    }

    @Override
    public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<CortexVertex, CortexEdge> g, int junctions, int size, int edges, int childColor, Set<Integer> parentColors) {
        return false;
    }

    boolean foundNovels = false;
    boolean arrivedOnShore = false;
    List<CortexVertex> novels = new ArrayList<>();
    int distanceFromLastNovel = 0;
    int distanceFromShore = 0;
    boolean reachedShore = false;

    @Override
    public boolean hasTraversalSucceeded(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedGraph<CortexVertex, CortexEdge> previousGraph, CortexGraph rois) {
        if (foundNovels) {
            distanceFromLastNovel++;
        }

        if (rois.findRecord(cv.getCr().getCortexKmer()) != null) {
            foundNovels = true;
            distanceFromLastNovel++;
            novels.add(cv);
            distanceFromShore = 0;
        } else {
            distanceFromShore++;
        }

        /*
        if (foundNovels) {
            for (int joiningColor : joiningColors) {
                if (cv.getCr().getCoverage(joiningColor) > 0) {
                    reachedShore = true;
                }
            }
        }

        return reachedShore;
        */

        return foundNovels && (distanceFromShore >= 10 || currentTraversalDepth >= 2 || numAdjacentEdges == 0 || childrenAlreadyTraversed);
    }

    @Override
    public boolean hasTraversalFailed(CortexVertex cv, boolean goForward, int traversalColor, Set<Integer> joiningColors, int currentTraversalDepth, int currentGraphSize, int numAdjacentEdges, boolean childrenAlreadyTraversed, DirectedGraph<CortexVertex, CortexEdge> previousGraph, CortexGraph rois) {
        return !foundNovels && (currentGraphSize >= 100 || currentTraversalDepth >= 2 || numAdjacentEdges == 0);
    }
}
