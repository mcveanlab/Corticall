package uk.ac.ox.well.indiana.analyses.roi;

import com.google.common.base.Joiner;
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
import java.util.*;

public class VisualizeSequenceTree extends Module {
    @Argument(fullName="sequences", shortName="s", doc="Sequences FASTA")
    public FastaSequenceFile FASTA;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Argument(fullName="skipSimplification", shortName="ss", doc="Don't simplify the graph")
    public Boolean SKIP_SIMPLIFICATION = false;

    @Argument(fullName="highlight", shortName="hl", doc="Colon-separated regex:color list (e.g. 'Pf3D7:#ff0000') specifying the sequence(s) to highlight with the specified color.", required=false)
    public ArrayList<String> HIGHLIGHT;

    @Argument(fullName="expansionFactor", shortName="ef", doc="Expansion factor for edge thickness")
    public Float EXPANSION_FACTOR = 0.1f;

    @Output
    public File out;

    public class CustomLabelProvider implements VertexNameProvider<String> {
        @Override
        public String getVertexName(String vertex) {
            return "";
            //return vertex;
            //return String.valueOf(vertex.length());
            //return SKIP_SIMPLIFICATION ? vertex : String.valueOf(vertex.length());
        }
    }

    public class CustomVertexAttributeProvider implements ComponentAttributeProvider<String> {
        @Override
        public Map<String, String> getComponentAttributes(String s) {
            Map<String, String> attrs = new HashMap<String, String>();
            attrs.put("shape", "circle");
            attrs.put("width", String.valueOf(Math.log10(s.length() - KMER_SIZE + 1)/10.0 + 0.12));
            attrs.put("fontsize", "1");

            return attrs;
        }
    }

    public class CustomEdgeAttributeProvider implements ComponentAttributeProvider<DefaultEdge> {
        private Map<DefaultEdge, Map<String, Float>> edgeWeights = new HashMap<DefaultEdge, Map<String, Float>>();
        private Map<DefaultEdge, Set<String>> edgeColors = new HashMap<DefaultEdge, Set<String>>();

        public void incrementWeight(DefaultEdge edge, String color) {
            if (!edgeWeights.containsKey(edge)) {
                edgeWeights.put(edge, new HashMap<String, Float>());
            }

            if (!edgeWeights.get(edge).containsKey(color)) {
                edgeWeights.get(edge).put(color, 1.0f);
            } else {
                edgeWeights.get(edge).put(color, edgeWeights.get(edge).get(color) + EXPANSION_FACTOR);
            }
        }

        public void addColor(DefaultEdge edge, String color) {
            if (!edgeColors.containsKey(edge)) {
                edgeColors.put(edge, new HashSet<String>());
            }

            edgeColors.get(edge).add(color);
        }

        @Override
        public Map<String, String> getComponentAttributes(DefaultEdge edge) {
            Map<String, String> attrs = new HashMap<String, String>();
            attrs.put("arrowhead", "none");

            Map<String, Integer> colorWidths = new HashMap<String, Integer>();

            if (edgeColors.containsKey(edge)) {
                attrs.put("color", Joiner.on(":").join(edgeColors.get(edge)));

                for (String color : edgeColors.get(edge)) {
                    if (!colorWidths.containsKey(color)) {
                        colorWidths.put(color, 1);
                    } else {
                        colorWidths.put(color, colorWidths.get(color));
                    }
                }

                List<Float> widths = new ArrayList<Float>();
                for (String color : edgeColors.get(edge)) {
                    widths.add(edgeWeights.get(edge).get(color));
                }

                attrs.put("penwidth", Joiner.on(":").join(widths));
            }

            return attrs;
        }
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
        // Load genes
        DirectedGraph<String, DefaultEdge> graph = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);
        Map<String, String> sequences = new HashMap<String, String>();

        ReferenceSequence rseq;
        while ((rseq = FASTA.nextSequence()) != null) {
            // Load each gene into its own graph object.
            DirectedGraph<String, DefaultEdge> sampleGraph = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);
            Set<String> prevKmers = new HashSet<String>();
            String prevKmer = null;

            String seq = new String(rseq.getBases());

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

            // Add all the graphs together.
            Graphs.addGraph(graph, sampleGraph);

            sequences.put(rseq.getName(), seq);
        }

        // Find all of the branch points in the graph.
        Set<String> branchVertices = new HashSet<String>();
        for (String vertex : graph.vertexSet()) {
            if (graph.inDegreeOf(vertex) != 1 || graph.outDegreeOf(vertex) != 1) {
                branchVertices.add(vertex);
            }
        }

        DirectedGraph<String, DefaultEdge> finalGraph;
        if (SKIP_SIMPLIFICATION) {
            finalGraph = graph;
        } else {
            // Build the complete graph again, but this time, break the sequences into
            // variable sized chunks based on where the branching kmers are found.
            finalGraph = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);
            for (String seqName : sequences.keySet()) {
                String seq = sequences.get(seqName);

                StringBuilder contig = new StringBuilder();
                String prevContig = null;

                for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                    String kmer = seq.substring(i, i + KMER_SIZE);

                    if (contig.length() == 0) {
                        contig.append(kmer);
                    } else {
                        contig.append(kmer.charAt(kmer.length() - 1));
                    }

                    if (branchVertices.contains(kmer) && i > 0) {
                        finalGraph.addVertex(contig.toString());

                        if (prevContig != null) {
                            finalGraph.addEdge(prevContig, contig.toString());
                        }

                        prevContig = contig.toString();
                        contig = new StringBuilder();
                    }
                }
            }
        }

        // Annotate the graph with some useful information.
        CustomEdgeAttributeProvider eap = new CustomEdgeAttributeProvider();

        Map<String, String> colorMap = new HashMap<String, String>();
        if (HIGHLIGHT != null) {
            for (String regexColor : HIGHLIGHT) {
                String[] pieces = regexColor.split(":");
                String regex = pieces[0];
                String color = pieces[1];

                colorMap.put(regex, color);
            }
        }

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
                    boolean found = false;
                    for (String regex : colorMap.keySet()) {
                        if (seqName.matches(regex)) {
                            eap.addColor(edge, colorMap.get(regex));
                            found = true;

                            eap.incrementWeight(edge, colorMap.get(regex));
                            //eap.incrementWeight(edge, colorMap.get(regex));
                            //eap.incrementWeight(edge, colorMap.get(regex));
                        }
                    }

                    if (!found) {
                        eap.addColor(edge, "#000000");
                        eap.incrementWeight(edge, "#000000");
                    }
                }
            }
        }

        log.info("Sequences: {}", sequences.size());
        log.info("Vertices: {}", graph.vertexSet().size());
        log.info("Edges: {}", graph.edgeSet().size());
        log.info("Simplified vertices: {}", finalGraph.vertexSet().size());
        log.info("Simplified branch vertices: {}", branchVertices.size());
        log.info("Simplified edges: {}", finalGraph.edgeSet().size());

        // Write the graph to disk.
        writeGraph(finalGraph, eap, out);
    }
}
