package uk.ac.ox.well.indiana.utils.traversal;

import org.jetbrains.annotations.NotNull;
import org.jgrapht.DirectedGraph;
import org.jgrapht.traverse.TopologicalOrderIterator;
import org.junit.BeforeClass;
import org.testng.annotations.Test;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.*;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.stoppingconditions.DustStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.LongWalkStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.OrphanStopper;

import java.io.File;
import java.util.*;

/**
 * Created by kiran on 10/05/2017.
 */
public class TraversalEngineTest {
    private File tempDir = new File("testdata");
    private File simpleGraph = new File(tempDir.getAbsolutePath() + "/_simplegraph.ctx");

    @BeforeClass
    public void writeTemporaryTestFiles() {
        CortexGraphWriter cgw = new CortexGraphWriter(simpleGraph);
        cgw.setHeader(constructCortexHeader());

        String seqMother = "ATGTACCGTGTTGGCACGAAGTTTTGCTGATATAGGAGATATCGTACGCGGTAAAGATCTGTATCTCGGTAATCCACAAGAAAGTACACAAAGAATAATATTAGAAAATAATTTGAAAGATATTTTCGCGAAAATACATAGTGACGTGATGTCAACGAGCGGGAGTAATGGGAGGGCGCTACAAAAACGCTACAAAGATACTGATAATTATTATGAATTGAGAGAAGATTGGTGGGCACTTAATAGAGACCAAGTATGGAAAGCTATCACATGCAATGCTGGGGGTGGTAATAGATATTTTCGACAAACATGTGGTTCAGGAGAATGGGCTAAAGACAAATGCCGGTGTAAGGACGACAAGGTCCCCACATATTTTGACTATGTGCCACAGTATCTTCGCTGGTTCGAGGAATGGGCCGAAGATTTTTGTAGATTAAGGAAACATAAATTAAAAGATGCTAAAAACAAATGTCGTGGAGATAGTGGTAACGATAGATATTGTGATCTTAATAGGTATGATTGCACACAAACTATTAGAGGAAATGAACATTTTGTTGAAAAGGATGATTGTAAAGGTTGTCAGTATTCGTGCGCTCATTTTGTGAACTGGATAGATAACCAAAAACTAGAATTTGAAAAACAAAAAGAAAAATATACAAAAGAAATTAAAAAAAAGCATCCAACAACCATAATAATAA";
        String seqFather = "ATGTACCGTGTTGGCACGAAGTTTTGCTGATATAGGAGATATCGTACGCGATAAAGATCTGTATCTCGGTAATCCACAAGAAAGTACACAAAGAATAATATTAGAAAATAATTTGAAAGATATTTACGCGAAAATACATAGTGACGTGATGTCAACGAGCGGGAGTAATGGGAGGGCGCTACAAAAACGCTACAAAGATACTGATAATTATTATGAATTGAGAGACGATTGGTGGGCACTTAATAGAGACCAAGTATGGAAAGCTATCACATGCATTGCTGGGGGTGGTAATAGATATTTTCGACAAACATGTGGATCAGGAGAAAGGGCTAAAGACAAATGCCGGTGTAAGGACGACAAGGTCCAACGATTACCAGAACCATGGCCACAGTATCTTCGCTGGTTCGAGGAATGGGCCGAAGATTTTTGTAGAAAAAGGAAACACAAATTAAAAGATGCTAAAAACAAATGTCGTGGAGATAGTGGTAACGATAGATATTGTGATCTTAATAGGTATGATTGCACACAAACTATTAGAGGAAATCAACATTTTGTTGAAAAGGATGATTGTAAAGGTTGTCAGTATTCGTGCGCTCATTTTGTGAACTGGATAGATAACCAAAAACTAGAATTTGAAAAACAAACCGAAAAATATACAAAAGAAATTAAAAAAAAGCATCCAACAACCATAATAATAA";
        String seqChild  = "ATGTACCGTGTTGGCACGAAGTTTTGCTGATATAGGAGATATCGTACGCGGTAAAGATCTGTATCTCGGTAATCCACAAGAAAGTACACAAAGAATAATATTAGAAAATAATTTGAAAGATATTTTCGCGAAAATACATAGTGACGTGATGTCAACGAGCGGGAGTAATGGGAGGGCGCTACAAAAACGCTACAAAGATACTGATAATTATTATGAATTGAGAGAAGATTGGTGGGCACTTAATAGAGACCAAGTATGGAAAGCTATCACATGCAATGCTGGGGGTGGTAATAGATATTTTCGACAAACATGTGGTTCAGGAGAAAGGGCTAAAGACAAATGCCGGTGTAAGGACGACAAGGTCCAACGATTACCAGAACCATGGCCACAGTATCTTCGCTGGTTCGAGGAATGGGCCGAAGATTTTTGTAGAAAAAGGAAACACAAATTAAAAGATGCTAAAAACAAATGTCGTGGAGATAGTGGTAACGATAGATATTGTGATCTTAATAGGTATGATTGCACACAAACTATTAGAGGAAATCAACATTTTGTTGAAAAGGATGATTGTAAAGGTTGTCAGTATTCGTGCGCTCATTTTGTGAACTGGATAGATAACCAAAAACTAGAATTTGAAAAACAAACCGAAAAATATACAAAAGAAATTAAAAAAAAGCATCCAACAACCATAATAATAA";
        // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------mom-||-dad------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        Map<CortexKmer, CortexRecord> crs = new TreeMap<>();

        for (int i = 0; i <= seqMother.length() - cgw.getHeader().getKmerSize(); i++) {
            String sk = seqMother.substring(i, i + cgw.getHeader().getKmerSize());
            CortexKmer ck = new CortexKmer(sk);

            //CortexRecord cr = new CortexRecord(new CortexBinaryKmer(sk.getBytes()), );

        }
    }

    @NotNull
    private CortexHeader constructCortexHeader() {
        CortexHeader ch = new CortexHeader();
        ch.setVersion(6);
        ch.setNumColors(3);
        ch.setKmerSize(31);
        ch.setKmerBits(CortexUtils.getKmerBits(31));

        CortexColor cmother = new CortexColor();
        cmother.setSampleName("mother");

        CortexColor cfather = new CortexColor();
        cfather.setSampleName("father");

        CortexColor cchild = new CortexColor();
        cchild.setSampleName("child");

        ch.addColor(cmother);
        ch.addColor(cfather);
        ch.addColor(cchild);

        return ch;
    }

    @Test
    public void testTraversalEngine() {
        CortexGraph g = new CortexGraph("3D7xHB3.k47.clean.infer.ctx");

        int childColor = g.getColorForSampleName("PG0051-C.ERR019061.md");
        Set<Integer> parentColors = new HashSet<>(g.getColorsForSampleNames(Arrays.asList("PG0051-C.ERR019061.md", "PG0052-C.ERR019054.md")));
        Set<Integer> refColors = new HashSet<>(g.getColorsForSampleNames(Arrays.asList("3D7", "HB3")));

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(childColor)
                .joiningColors(parentColors)
                .recruitmentColors(refColors)
                .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.AND)
                .traversalDirection(TraversalEngineConfiguration.TraversalDirection.BOTH)
                .connectUnusedNeighbors(false)
                .stopper(new LongWalkStopper())
                .graph(g)
                .make();

        CortexKmer ck = new CortexKmer("ATGGAAATTGTATAAATAAAGATAATGACAACACATGTAAAAATTCA");
        String sk = ck.getKmerAsString();

        DirectedGraph<CortexVertex, CortexEdge> walkNew = e.dfs(sk);

        TopologicalOrderIterator<CortexVertex, CortexEdge> toi = new TopologicalOrderIterator<>(walkNew);

        StringBuilder sb = new StringBuilder();
        while (toi.hasNext()) {
            CortexVertex cv = toi.next();
            if (sb.length() == 0) {
                sb.append(cv.getSk());
            } else {
                sb.append(cv.getSk().substring(cv.getSk().length() - 1));
            }
        }

        System.out.println(sb);
        System.out.println();
    }
}
