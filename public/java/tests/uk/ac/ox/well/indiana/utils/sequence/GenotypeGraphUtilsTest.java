package uk.ac.ox.well.indiana.utils.sequence;

import org.apache.commons.math3.util.Pair;
import org.testng.annotations.Test;

public class GenotypeGraphUtilsTest {
    //DirectedGraph<AnnotatedVertex, AnnotatedEdge> ag = GenotypeGraphUtils.loadLocalSubgraph(stretch, CLEAN, GRAPHS, novelKmers);

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

        return new Pair<>(ref, alt);
    }

    @Test(enabled = false)
    public void testSNPRecovery() {
        Pair<String, String> h = generateSNPHaplotypes(10000);
    }
}
