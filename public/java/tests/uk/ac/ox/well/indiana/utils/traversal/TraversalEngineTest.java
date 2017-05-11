package uk.ac.ox.well.indiana.utils.traversal;

import org.jgrapht.DirectedGraph;
import org.jgrapht.traverse.TopologicalOrderIterator;
import org.testng.annotations.Test;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.stoppingconditions.DustStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.LongWalkStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.OrphanStopper;

import java.util.*;

/**
 * Created by kiran on 10/05/2017.
 */
public class TraversalEngineTest {
    @Test
    public void testTraversalEngine() {
        CortexGraph g = new CortexGraph("3D7xHB3.k47.clean.infer.ctx");

        int childColor = g.getColorForSampleName("PG0051-C.ERR019061.md");
        Set<Integer> parentColors = new HashSet<>(g.getColorsForSampleNames(Arrays.asList("PG0051-C.ERR019061.md", "PG0052-C.ERR019054.md")));
        Set<Integer> refColors = new HashSet<>(g.getColorsForSampleNames(Arrays.asList("3D7", "HB3")));
        Set<Integer> ignoreColors = new HashSet<>(g.getColorsForSampleNames(Arrays.asList("PG0445-C", "PG0453-C", "PG0457-C", "PG0477-C")));

        List<Integer> displayColors = new ArrayList<>();
        displayColors.add(childColor);
        displayColors.addAll(parentColors);
        displayColors.addAll(refColors);

        Set<CortexKmer> novelKmers = new HashSet<>();
        novelKmers.add(new CortexKmer("ATGGAAATTGTATAAATAAAGATAATGACAACACATGTAAAAATTCA"));

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(childColor)
                .joiningColors(parentColors)
                .recruitmentColors(refColors)
                .displayColors(displayColors)
                .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.AND)
                .traversalDirection(TraversalEngineConfiguration.TraversalDirection.BOTH)
                .connectUnusedNeighbors(false)
                .stopper(new LongWalkStopper())
                .graph(g)
                .make();

        for (CortexKmer ck : novelKmers) {
            String kmer = ck.getKmerAsString();

            DirectedGraph<CortexVertex, CortexEdge> walkNew = e.dfs(kmer);

            TopologicalOrderIterator<CortexVertex, CortexEdge> toi = new TopologicalOrderIterator<>(walkNew);

            while (toi.hasNext()) {
                System.out.println(toi.next());
            }

            System.out.println();
        }
    }
}
