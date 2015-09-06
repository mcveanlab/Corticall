package uk.ac.ox.well.indiana.commands.gg;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

public class ParentTraversalStopper extends AbstractTraversalStopper<AnnotatedVertex, AnnotatedEdge> {
    private boolean sawPredecessorFirst = false;
    private boolean sawSuccessorFirst = false;

    private boolean isLowComplexity(CortexRecord cr) {
        byte edges[][] = cr.getEdgesAsBytes();

        int numEdges = 0;
        for (int e = 0; e < 8; e++) {
            if (edges[0][e] != '.') {
                numEdges++;
            }
        }

        //return numEdges > 6;

        return false;
    }

    @Override
    public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions, int size, int edges) {
        String fw = cr.getKmerAsString();
        String rc = SequenceUtils.reverseComplement(fw);

        AnnotatedVertex v = null;

        boolean rejoinedGraph = false;

        if (g.containsVertex(new AnnotatedVertex(fw)) || g.containsVertex(new AnnotatedVertex(rc)) || g.containsVertex(new AnnotatedVertex(fw, true)) || g.containsVertex(new AnnotatedVertex(rc, true))) {
            rejoinedGraph = true;

            for (AnnotatedVertex av : g.vertexSet()) {
                /*
                if ((av.getKmer().equals(fw) || av.getKmer().equals(rc)) && (av.flagIsSet("predecessor") || av.flagIsSet("successor"))) {
                    rejoinedGraph = false;

                    break;
                }
                */

                if (av.getKmer().equals(fw) || av.getKmer().equals(rc)) {
                    if (!sawPredecessorFirst && !sawSuccessorFirst) {
                        if (av.flagIsSet("predecessor")) { sawPredecessorFirst = true; }
                        else if (av.flagIsSet("successor")) { sawSuccessorFirst = true; }
                    }

                    if ((sawPredecessorFirst && av.flagIsSet("predecessor")) || (sawSuccessorFirst && av.flagIsSet("successor"))) {
                        rejoinedGraph = false;
                        break;
                    }
                }
            }
        }

        //return size > 1 && (g.containsVertex(new AnnotatedVertex(fw)) || g.containsVertex(new AnnotatedVertex(rc)) || g.containsVertex(new AnnotatedVertex(fw, true)) || g.containsVertex(new AnnotatedVertex(rc, true)));
        return size > 1 && rejoinedGraph;
    }

    @Override
    public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions, int size, int edges) {
        return size > 1000 || junctions >= maxJunctionsAllowed() || isLowComplexity(cr);
    }

    @Override
    public int maxJunctionsAllowed() {
        return 5;
    }
}
