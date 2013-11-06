package uk.ac.ox.well.indiana.analyses.roi;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.ext.DOTExporter;
import org.jgrapht.ext.StringNameProvider;
import org.jgrapht.ext.VertexNameProvider;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.tools.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

public class VisualizeSequenceTree extends Module {
    @Argument(fullName="sequences", shortName="s", doc="Sequences FASTA")
    public FastaSequenceFile FASTA;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Output
    public File out;

    public class CustomLabelProvider implements VertexNameProvider<String> {
        @Override
        public String getVertexName(String vertex) {
            return vertex.length() > KMER_SIZE ? String.valueOf(vertex.length()) : vertex;
        }
    }

    private Set<String> firstVertices(DirectedGraph<String, DefaultEdge> graph) {
        Set<String> zeroInDegrees = new HashSet<String>();

        for (String vertex : graph.vertexSet()) {
            if (graph.inDegreeOf(vertex) == 0) {
                zeroInDegrees.add(vertex);
            }
        }

        return zeroInDegrees;
    }

    private DirectedGraph<String, DefaultEdge> simplify(DirectedGraph<String, DefaultEdge> graph, String startingVertex) {
        DirectedGraph<String, DefaultEdge> sg1 = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);

        StringBuilder contig = new StringBuilder(startingVertex.substring(0, startingVertex.length() - 1));

        Set<String> seenVertices = new HashSet<String>();
        String currentVertex = startingVertex;

        while (graph.outDegreeOf(currentVertex) > 0 && !seenVertices.contains(currentVertex)) {
            contig.append(currentVertex.charAt(currentVertex.length() - 1));

            seenVertices.add(currentVertex);

            Set<DefaultEdge> edges = graph.outgoingEdgesOf(currentVertex);

            if (edges.size() == 1) {
                DefaultEdge edge = edges.iterator().next();

                String nextVertex = graph.getEdgeTarget(edge);

                if (graph.inDegreeOf(nextVertex) == 1) {
                    currentVertex = nextVertex;
                } else {
                    Graphs.addGraph(sg1, graph);
                    sg1.removeAllVertices(seenVertices);
                    sg1.addVertex(contig.toString());

                    if (graph.inDegreeOf(startingVertex) > 0) {
                        for (DefaultEdge inEdge : graph.incomingEdgesOf(startingVertex)) {
                            sg1.addEdge(graph.getEdgeSource(inEdge), contig.toString());
                        }
                    }

                    sg1.addEdge(contig.toString(), nextVertex);

                    return simplify(sg1, nextVertex);
                }
            } else {
                Graphs.addGraph(sg1, graph);
                sg1.removeAllVertices(seenVertices);
                sg1.addVertex(contig.toString());

                if (graph.inDegreeOf(startingVertex) > 0) {
                    for (DefaultEdge inEdge : graph.incomingEdgesOf(startingVertex)) {
                        sg1.addEdge(graph.getEdgeSource(inEdge), contig.toString());
                    }
                }

                for (DefaultEdge edge : edges) {
                    String nextVertex = graph.getEdgeTarget(edge);

                    sg1.addEdge(contig.toString(), nextVertex);

                    sg1 = simplify(sg1, nextVertex);
                }

                return sg1;
            }
        }

        seenVertices.add(currentVertex);

        contig.append(currentVertex.charAt(currentVertex.length() - 1));

        Graphs.addGraph(sg1, graph);
        sg1.removeAllVertices(seenVertices);
        sg1.addVertex(contig.toString());

        if (graph.inDegreeOf(startingVertex) > 0) {
            for (DefaultEdge inEdge : graph.incomingEdgesOf(startingVertex)) {
                sg1.addEdge(graph.getEdgeSource(inEdge), contig.toString());
            }
        }

        return sg1;
    }

    private DirectedGraph<String, DefaultEdge> simplify(DirectedGraph<String, DefaultEdge> graph) {
        Set<String> firstVertices = firstVertices(graph);

        DirectedGraph<String, DefaultEdge> simplifiedGraph = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);
        Graphs.addGraph(simplifiedGraph, graph);

        for (String vertex : firstVertices) {
            simplifiedGraph = simplify(simplifiedGraph, vertex);
        }

        return simplifiedGraph;
    }

    public void writeGraph(DirectedGraph<String, DefaultEdge> g, File f) {
        DOTExporter<String, DefaultEdge> exporter = new DOTExporter<String, DefaultEdge>(new StringNameProvider<String>(), new CustomLabelProvider(), null, null, null);

        try {
            exporter.export(new FileWriter(f), g);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void execute() {
        DirectedGraph<String, DefaultEdge> graph = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);

        String startKmer = null;

        ReferenceSequence rseq;
        while ((rseq = FASTA.nextSequence()) != null) {
            String seq = new String(rseq.getBases());

            String prevKmer = null;
            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                String kmer = seq.substring(i, i + KMER_SIZE);

                if (startKmer == null) {
                    startKmer = kmer;
                }

                graph.addVertex(kmer);

                if (prevKmer != null) {
                    graph.addEdge(prevKmer, kmer);
                }

                prevKmer = kmer;
            }
        }

        DirectedGraph<String, DefaultEdge> simplifiedGraph = simplify(graph);

        writeGraph(simplifiedGraph, out);
    }
}
