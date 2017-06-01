package uk.ac.ox.well.indiana.utils.visualizer;

import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.testng.annotations.Test;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;

public class GraphVisualizerTest {
    @Test
    public void testVisualizer() {
        DirectedGraph<CortexVertex, CortexEdge> g = new DefaultDirectedGraph<>(CortexEdge.class);

        String template = new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(100));
        StringBuilder sb1 = new StringBuilder(template);
        sb1.setCharAt(50, 'A');
        StringBuilder sb2 = new StringBuilder(template);
        sb2.setCharAt(50, 'C');

        String seq1 = sb1.toString();
        String seq2 = sb2.toString();

        String[] seqs = new String[] { seq1, seq2 };

        int kmerSize = 47;

        for (int q = 0; q < seqs.length; q++) {
            String seq = seqs[q];

            CortexVertex lv = null;

            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                String kmer = seq.substring(i, i + kmerSize);

                CortexVertex cv = new CortexVertex(kmer, null, 0);
                g.addVertex(cv);

                if (lv != null) {
                    g.addEdge(lv, cv, new CortexEdge(q, 1.0));
                }

                lv = cv;
            }
        }

        GraphVisualizer gv = new GraphVisualizer(9000);

        //int response = gv.display(g);

        //System.out.println("Response: " + response);
    }
}
