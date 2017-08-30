package uk.ac.ox.well.cortexjdk.utils.caller;

import org.jgrapht.GraphPath;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.DestinationStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.*;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Created by kiran on 30/07/2017.
 */
public class BubbleCaller {
    private BubbleCaller() {}

    static public Bubble call(TraversalEngineConfiguration config, List<CortexVertex> lwalk, int refColor, int altColor, CortexKmer novelKmer) {
        TraversalEngine e = new TraversalEngineFactory()
            .configuration(config)
            .traversalColor(refColor)
            .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.OR)
            .traversalDirection(TraversalEngineConfiguration.TraversalDirection.BOTH)
            .stoppingRule(DestinationStopper.class)
            .make();

        DirectedWeightedPseudograph<CortexVertex, CortexEdge> sg = TraversalEngine.toGraph(config.getGraph(), lwalk, altColor, refColor);

        Set<CortexVertex> iss = new HashSet<>();
        Set<CortexVertex> oss = new HashSet<>();

        for (CortexVertex v : lwalk) {
            if (TraversalEngine.inDegree(sg, v) > 1) {
                iss.add(v);
            }

            if (TraversalEngine.outDegree(sg, v) > 1) {
                oss.add(v);
            }
        }

        if (iss.size() == 0 || oss.size() == 0) {
            return null;
        }

        for (CortexVertex is : iss) {
            for (CortexVertex os : oss) {
                DirectedWeightedPseudograph<CortexVertex, CortexEdge> p = new DirectedWeightedPseudograph<>(CortexEdge.class);
                p.addVertex(os);

                e.getConfiguration().setPreviousTraversal(p);

                DirectedWeightedPseudograph<CortexVertex, CortexEdge> w = e.dfs(is.getSk());

                if (w != null) {
                    Graphs.addGraph(sg, w);
                }
            }
        }

        PathFinder dspRef = new PathFinder(sg, refColor);
        PathFinder dspAlt = new PathFinder(sg, altColor);

        Set<Bubble> bubbles = new HashSet<>();

        for (CortexVertex os : oss) {
            for (CortexVertex is : iss) {
                if (!os.equals(is)) {
                    GraphPath<CortexVertex, CortexEdge> pRef = dspRef.getPathFinder(os, is);
                    GraphPath<CortexVertex, CortexEdge> pAlt = dspAlt.getPathFinder(os, is, novelKmer, true);

                    bubbles.add(new Bubble(pRef, pAlt, null, null));
                }
            }
        }

        return bubbles.iterator().next();
    }
}
