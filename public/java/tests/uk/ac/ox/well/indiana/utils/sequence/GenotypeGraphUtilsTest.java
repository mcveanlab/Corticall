package uk.ac.ox.well.indiana.utils.sequence;

import com.google.common.base.Joiner;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.DirectedGraph;
import org.jgrapht.GraphPath;
import org.jgrapht.alg.DijkstraShortestPath;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.testng.annotations.Test;
import uk.ac.ox.well.indiana.commands.gg.AnnotatedEdge;
import uk.ac.ox.well.indiana.commands.gg.AnnotatedVertex;
import uk.ac.ox.well.indiana.commands.gg.GenotypeGraphUtils;

import java.util.*;

public class GenotypeGraphUtilsTest {
    //DirectedGraph<AnnotatedVertex, AnnotatedEdge> ag = GenotypeGraphUtils.loadLocalSubgraph(stretch, CLEAN, DIRTY, novelKmers);

    private Pair<String, String> generateSNPHaplotypes(int length) {
        int mid = length / 2;
        String ref = new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(length));

        StringBuilder altb = new StringBuilder(ref);
        char newbase = ref.charAt(mid);
        while (newbase == ref.charAt(mid)) {
            newbase = (char) SequenceUtils.generateRandomNucleotideSequenceOfLengthN(1)[0];
        }

        altb.setCharAt(mid, newbase);

        String alt = altb.toString();

        return new Pair<String, String>(ref, alt);
    }

    @Test(enabled = false)
    public void testSNPRecovery() {
        Pair<String, String> h = generateSNPHaplotypes(10000);
    }
}
