package uk.ac.ox.well.indiana.commands.gg;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

public class ParentTraversalStopper extends AbstractTraversalStopper<AnnotatedVertex, AnnotatedEdge> {
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
    @Override
    public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions, int size, int edges) {
        String fw = cr.getKmerAsString();
        String rc = SequenceUtils.reverseComplement(fw);

        return size > 1 && (g.containsVertex(new AnnotatedVertex(fw)) || g.containsVertex(new AnnotatedVertex(rc)) || g.containsVertex(new AnnotatedVertex(fw, true)) || g.containsVertex(new AnnotatedVertex(rc, true)));
    }

    @Override
    public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions, int size, int edges) {
        return size > 100 || junctions >= maxJunctionsAllowed() || isLowComplexity(cr);
    }

    @Override
    public int maxJunctionsAllowed() {
        return 3;
    }
}
