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
import uk.ac.ox.well.indiana.utils.io.jgrapht.MultiWeightEdge;

import java.io.*;
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

    @Argument(fullName="minThickness", shortName="minThickness", doc="Minimum thickness of plotted elements")
    public Float MIN_THICKNESS = 0.3f;

    @Argument(fullName="maxThickness", shortName="maxThickness", doc="Maximum thickness of plotted elements")
    public Float MAX_THICKNESS = 20.0f;

    @Argument(fullName="withWeight", shortName="ww", doc="Enable weights")
    public Boolean WITH_WEIGHT = false;

    @Argument(fullName="withVertexLabels", shortName="wvl", doc="Enable vertex labels")
    public Boolean WITH_VERTEX_LABELS = false;

    @Argument(fullName="withEdgeLabels", shortName="wel", doc="Enable edge labels")
    public Boolean WITH_EDGE_LABELS = false;

    @Output
    public File out;

    private int numSequences = 0;

    public class CustomLabelProvider implements VertexNameProvider<String> {
        @Override
        public String getVertexName(String vertex) {
            //return "";
            return vertex;
            //return String.valueOf(vertex.length());
            //return SKIP_SIMPLIFICATION ? vertex : String.valueOf(vertex.length());
        }
    }

    public class CustomVertexAttributeProvider implements ComponentAttributeProvider<String> {
        private Set<String> branchVertices;

        public CustomVertexAttributeProvider(Set<String> branchVertices) {
            this.branchVertices = branchVertices;
        }

        @Override
        public Map<String, String> getComponentAttributes(String s) {
            Map<String, String> attrs = new HashMap<String, String>();
            attrs.put("shape", "circle");
            attrs.put("width", String.valueOf(Math.log10(s.length() - KMER_SIZE + 1)/10.0 + 0.12));

            //if (SKIP_SIMPLIFICATION) {
                attrs.put("color", branchVertices.contains(s) ? "red" : "black");
            //}

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
                edgeWeights.get(edge).put(color, edgeWeights.get(edge).get(color) + 1.0f);
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
            attrs.put("arrowhead", "normal");
            attrs.put("arrowsize", "0.8");

            if (edgeColors.containsKey(edge)) {
                attrs.put("color", Joiner.on(":").join(edgeColors.get(edge)));

                List<String> labels = new ArrayList<String>();
                List<Float> widths = new ArrayList<Float>();
                for (String color : edgeColors.get(edge)) {
                    float weight = edgeWeights.get(edge).get(color);
                    float scaledWeight = MIN_THICKNESS + weight*((MAX_THICKNESS - MIN_THICKNESS) / (numSequences - 1));

                    widths.add(scaledWeight);
                    labels.add(String.format("%d", (int) weight));
                }

                attrs.put("penwidth", Joiner.on(":").join(widths));
                attrs.put("label", Joiner.on(":").join(labels));
                attrs.put("weight", Joiner.on(":").join(labels));
            }

            return attrs;
        }
    }

    public void writeGraph(DirectedGraph<String, MultiWeightEdge> g, File f) {
        try {
            PrintStream ps = new PrintStream(f);

            String indent = "  ";

            ps.println("digraph G {");
            ps.println(indent + "rankdir=\"LR\";");

            for (String vertex : g.vertexSet()) {
                if (WITH_VERTEX_LABELS) {
                    ps.println(indent + vertex + " [ label=\"" + vertex + "\" color=\"black\" width=\"0.10\" shape=\"circle\" fontsize=\"1\" ];");
                } else {
                    ps.println(indent + vertex + " [ label=\"\" color=\"black\" width=\"0.10\" shape=\"circle\" fontsize=\"1\" ];");
                }
                //ps.println(indent + vertex + " [ label=" + vertex + " color=\"black\" height=\"0.40\" shape=\"rect\" fontsize=\"1\" ];");
            }

            for (MultiWeightEdge e : g.edgeSet()) {
                Map<String, Integer> weights = e.getWeights();

                for (String color : weights.keySet()) {
                    float penwidth = MIN_THICKNESS + weights.get(color)*((MAX_THICKNESS - MIN_THICKNESS) / (numSequences - 1));

                    if (WITH_WEIGHT) {
                        if (WITH_EDGE_LABELS) {
                            ps.println(indent + g.getEdgeSource(e) + " -> " + g.getEdgeTarget(e) + " [ color=\"" + color + "\" label=\"" + weights.get(color) + "\" weight=\"" + weights.get(color) + "\" penwidth=\"" + penwidth + "\" arrowsize=\"0.8\" arrowhead=\"normal\" ];");
                        } else {
                            ps.println(indent + g.getEdgeSource(e) + " -> " + g.getEdgeTarget(e) + " [ color=\"" + color + "\" weight=\"" + weights.get(color) + "\" penwidth=\"" + penwidth + "\" arrowsize=\"0.8\" arrowhead=\"normal\" ];");
                        }
                    } else {
                        if (WITH_EDGE_LABELS) {
                            ps.println(indent + g.getEdgeSource(e) + " -> " + g.getEdgeTarget(e) + " [ color=\"" + color + "\" label=\"" + weights.get(color) + "\" penwidth=\"" + penwidth + "\" arrowsize=\"0.8\" arrowhead=\"normal\" ];");
                        } else {
                            ps.println(indent + g.getEdgeSource(e) + " -> " + g.getEdgeTarget(e) + " [ color=\"" + color + "\" penwidth=\"" + penwidth + "\" arrowsize=\"0.8\" arrowhead=\"normal\" ];");
                        }
                    }
                }
            }

            ps.println("}");
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
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
            numSequences++;
        }

        // Populate regex to color map
        Map<String, String> regexColorMap = new HashMap<String, String>();
        if (HIGHLIGHT != null) {
            for (String regexColor : HIGHLIGHT) {
                String[] pieces = regexColor.split(":");
                String regex = pieces[0];
                String color = pieces[1];

                regexColorMap.put(regex, color);
            }
        }

        // Populate seq to color map
        Map<String, String> seqColorMap = new HashMap<String, String>();
        for (String seqName : sequences.keySet()) {
            String color = "#000000";
            if (HIGHLIGHT != null) {
                for (String regex : regexColorMap.keySet()) {
                    if (seqName.matches(regex) || seqName.contains(regex)) {
                        color = regexColorMap.get(regex);
                    }
                }
            }

            seqColorMap.put(seqName, color);
        }

        // Find all of the branch points in the graph.
        Set<String> branchVertices = new HashSet<String>();
        for (String vertex : graph.vertexSet()) {
            if (graph.inDegreeOf(vertex) != 1 || graph.outDegreeOf(vertex) != 1) {
                branchVertices.add(vertex);
            }
        }

        DirectedGraph<String, MultiWeightEdge> finalGraph;
        if (SKIP_SIMPLIFICATION) {
            //finalGraph = graph;
            finalGraph = null;
        } else {
            // Build the complete graph again, but this time, break the sequences into
            // variable sized chunks based on where the branching kmers are found.
            finalGraph = new DefaultDirectedGraph<String, MultiWeightEdge>(MultiWeightEdge.class);

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
                            finalGraph.getEdge(prevContig, contig.toString()).incrementWeight(seqColorMap.get(seqName));
                        }

                        prevContig = contig.toString();
                        contig = new StringBuilder();
                        contig.append(kmer);
                    }
                }
            }
        }

        log.info("Sequences: {}", sequences.size());
        log.info("Edges: {}", graph.edgeSet().size());
        log.info("Vertices: {}", graph.vertexSet().size());
        log.info("Branch vertices: {}", branchVertices.size());
        log.info("Simplified vertices: {}", finalGraph.vertexSet().size());
        log.info("Simplified edges: {}", finalGraph.edgeSet().size());

        // Write the graph to disk.
        writeGraph(finalGraph, out);
    }
}
