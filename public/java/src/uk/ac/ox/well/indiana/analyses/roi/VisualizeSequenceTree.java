package uk.ac.ox.well.indiana.analyses.roi;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.ext.ComponentAttributeProvider;
import org.jgrapht.ext.StringNameProvider;
import org.jgrapht.ext.VertexNameProvider;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.tools.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.jgrapht.ExtendedDOTExporter;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class VisualizeSequenceTree extends Module {
    @Argument(fullName="sequences", shortName="s", doc="Sequences FASTA")
    public FastaSequenceFile FASTA;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Argument(fullName="skipSimplification", shortName="ss", doc="Don't simplify the graph")
    public Boolean SKIP_SIMPLIFICATION = false;

    @Output
    public File out;

    public class CustomLabelProvider implements VertexNameProvider<String> {
        @Override
        public String getVertexName(String vertex) {
            return "";
            //return String.valueOf(vertex.length());
            //return SKIP_SIMPLIFICATION ? vertex : String.valueOf(vertex.length());
        }
    }

    public class CustomVertexAttributeProvider implements ComponentAttributeProvider<String> {
        @Override
        public Map<String, String> getComponentAttributes(String s) {
            Map<String, String> attrs = new HashMap<String, String>();
            attrs.put("shape", "circle");
            attrs.put("height", "0.12");
            attrs.put("width", "0.12");
            attrs.put("fontsize", "1");

            return attrs;
        }
    }

    public class CustomEdgeAttributeProvider implements ComponentAttributeProvider<DefaultEdge> {
        private Map<DefaultEdge, Integer> edgeWeights = new HashMap<DefaultEdge, Integer>();

        public void incrementWeight(DefaultEdge edge) {
            if (edgeWeights.containsKey(edge)) {
                edgeWeights.put(edge, edgeWeights.get(edge) + 2);
            } else {
                edgeWeights.put(edge, 1);
            }
        }

        public int getWeight(DefaultEdge edge) {
            return edgeWeights.containsKey(edge) ? edgeWeights.get(edge) : 0;
        }

        @Override
        public Map<String, String> getComponentAttributes(DefaultEdge edge) {
            Map<String, String> attrs = new HashMap<String, String>();
            attrs.put("arrowhead", "none");

            if (edgeWeights.containsKey(edge)) {
                attrs.put("penwidth", String.valueOf(edgeWeights.get(edge)));
            }

            return attrs;
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
                        //if (sg1.containsVertex(graph.getEdgeSource(inEdge)) && sg1.containsVertex(contig.toString())) {
                            sg1.addEdge(graph.getEdgeSource(inEdge), contig.toString());
                        //}
                    }
                }

                for (DefaultEdge edge : edges) {
                    String nextVertex = graph.getEdgeTarget(edge);

                    //if (sg1.containsVertex(contig.toString()) && sg1.containsVertex(nextVertex)) {
                        sg1.addEdge(contig.toString(), nextVertex);
                    //}

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

        int count = 1;
        for (String vertex : firstVertices) {
            log.info("Simplifying graph ({}/{})", count, firstVertices.size());
            simplifiedGraph = simplify(simplifiedGraph, vertex);

            count++;
        }

        return simplifiedGraph;
    }

    public void writeGraph(DirectedGraph<String, DefaultEdge> g, CustomEdgeAttributeProvider eap, File f) {
        Map<String, String> attributes = new HashMap<String, String>();
        attributes.put("rankdir", "LR");
        ExtendedDOTExporter<String, DefaultEdge> exporter = new ExtendedDOTExporter<String, DefaultEdge>(
                new StringNameProvider<String>(),
                new CustomLabelProvider(),
                null,
                new CustomVertexAttributeProvider(),
                eap,
                attributes
        );

        try {
            exporter.export(new FileWriter(f), g);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void execute() {
        DirectedGraph<String, DefaultEdge> graph = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);
        Map<String, String> sequences = new HashMap<String, String>();

        ReferenceSequence rseq;
        while ((rseq = FASTA.nextSequence()) != null) {
            DirectedGraph<String, DefaultEdge> sampleGraph = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);

            String seq = new String(rseq.getBases());

            Set<String> prevKmers = new HashSet<String>();
            String prevKmer = null;

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                String kmer = seq.substring(i, i + KMER_SIZE);

                sampleGraph.addVertex(kmer);

                if (prevKmer != null && !prevKmers.contains(kmer)) {
                    sampleGraph.addEdge(prevKmer, kmer);
                }

                prevKmer = kmer;

                if (!SKIP_SIMPLIFICATION) {
                    prevKmers.add(kmer);
                }
            }

            Graphs.addGraph(graph, sampleGraph);

            sequences.put(rseq.getName(), seq);
        }

        DirectedGraph<String, DefaultEdge> finalGraph = SKIP_SIMPLIFICATION ? graph : simplify(graph);

        CustomEdgeAttributeProvider eap = new CustomEdgeAttributeProvider();

        for (String seqName : sequences.keySet()) {
            String seq = sequences.get(seqName);

            Set<String> vertices = new HashSet<String>();
            for (String vertex : finalGraph.vertexSet()) {
                if (seq.contains(vertex)) {
                    vertices.add(vertex);
                }
            }

            for (DefaultEdge edge : finalGraph.edgeSet()) {
                String sourceVertex = finalGraph.getEdgeSource(edge);
                String targetVertex = finalGraph.getEdgeTarget(edge);

                if (vertices.contains(sourceVertex) && vertices.contains(targetVertex)) {
                    int oldWeight = eap.getWeight(edge);

                    eap.incrementWeight(edge);

                    int newWeight = eap.getWeight(edge);
                }
            }
        }

        writeGraph(finalGraph, eap, out);
    }
}
