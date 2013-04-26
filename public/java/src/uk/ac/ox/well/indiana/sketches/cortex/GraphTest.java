package uk.ac.ox.well.indiana.sketches.cortex;

import org.jgrapht.DirectedGraph;
import org.jgrapht.ext.DOTExporter;
import org.jgrapht.ext.StringNameProvider;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.FileWriter;
import java.io.IOException;

public class GraphTest extends Tool {
    @Argument(fullName="g1", shortName="g1", doc="Graph 1")
    public String G1;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size", required=false)
    public Integer KMER_SIZE = 5;

    @Argument(fullName="canonicalize", shortName="c", doc="Canonicalize kmer?")
    public Boolean CANONICALIZE = false;

    @Override
    public void execute() {
        String seq = "TAGCTGTCTCTATGCTTTCTCTCTGCTCTATATATATAAAAAAATCGTTTCTGA";

        DirectedGraph<String, DefaultEdge> directedGraph = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);

        for (int i = 0; i < seq.length() - KMER_SIZE - 1; i++) {
            String kmer = seq.substring(i, i + KMER_SIZE);
            String nextKmer = seq.substring(i+1, i+1 + KMER_SIZE);

            if (CANONICALIZE) {
                kmer = new String(SequenceUtils.alphanumericallyLowestOrientation(kmer.getBytes()));
                nextKmer = new String(SequenceUtils.alphanumericallyLowestOrientation(nextKmer.getBytes()));
            }

            directedGraph.addVertex(kmer);
            directedGraph.addVertex(nextKmer);
            directedGraph.addEdge(kmer, nextKmer);
        }

        DOTExporter<String, DefaultEdge> exporter = new DOTExporter<String, DefaultEdge>(new StringNameProvider<String>(), null, null);

        try {
            exporter.export(new FileWriter(G1), directedGraph);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }
}
