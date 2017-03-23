package uk.ac.ox.well.indiana.utils.traversal;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.util.Set;

public class ChildTraversalStopper extends AbstractTraversalStopper<AnnotatedVertex, AnnotatedEdge> {
    private int goalSize = 0;
    private int goalDepth = 0;

    private boolean isNovel(CortexRecord cr, Set<Integer> childColors, Set<Integer> parentColors) {
        boolean isInParents = isSharedWithAParent(cr, parentColors);

        if (!isInParents) {
            for (int c : childColors) {
                if (cr.getCoverage(c) > 0) {
                    return true;
                }
            }
        }

        return false;
    }

    private boolean isSharedWithAParent(CortexRecord cr, Set<Integer> parentColors) {
        for (int c : parentColors) {
            if (cr.getCoverage(c) > 0) {
                return true;
            }
        }

        return false;
    }

    public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int depth, int size, int edges, Set<Integer> childColors, Set<Integer> parentColors) {
        if (goalSize == 0 && isSharedWithAParent(cr, parentColors)) {
            goalSize = size;
            goalDepth = depth;
        }

        if (goalSize > 0 && isNovel(cr, childColors, parentColors)) {
            goalSize = size;
            goalDepth = depth;
        }

        return (goalSize > 0 && (size >= goalSize + 50 || edges == 0));
    }

    @Override
    public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int depth, int size, int edges, Set<Integer> childColors, Set<Integer> parentColors) {
        return !isNovel(cr, childColors, parentColors) && (edges == 0 || depth >= maxJunctionsAllowed() || size > 5000);
    }

    @Override
    public int maxJunctionsAllowed() {
        return 5;
    }
}
