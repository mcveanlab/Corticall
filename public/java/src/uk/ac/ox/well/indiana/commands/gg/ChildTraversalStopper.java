package uk.ac.ox.well.indiana.commands.gg;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

public class ChildTraversalStopper extends AbstractTraversalStopper<AnnotatedVertex, AnnotatedEdge> {
    private int goalSize = 0;
    private int goalDepth = 0;

    private boolean isLowComplexity(CortexRecord cr) {
        byte edges[][] = cr.getEdgesAsBytes();

        int numEdges = 0;
        for (int e = 0; e < 8; e++) {
            if (edges[0][e] != '.') {
                numEdges++;
            }
        }

        return numEdges > 6;
    }

    public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int depth, int size, int edges) {
        if (goalSize == 0 && (cr.getCoverage(1) > 0 || cr.getCoverage(2) > 0)) {
            goalSize = size;
            goalDepth = depth;
        }

        if (goalSize > 0 && (cr.getCoverage(0) > 10 && cr.getCoverage(1) == 0 && cr.getCoverage(2) == 0)) {
            goalSize = size;
            goalDepth = depth;
        }

        return (goalSize > 0 && (size >= goalSize + 50 || isLowComplexity(cr) || edges == 0));
    }

    @Override
    public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int depth, int size, int edges) {
        //return (goalSize > 0 && depth >= goalDepth + 1) || edges == 0;
        return isLowComplexity(cr) || edges == 0 || depth >= maxJunctionsAllowed();
    }

    @Override
    public int maxJunctionsAllowed() {
        return 2;
    }
}
