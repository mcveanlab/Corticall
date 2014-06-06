package uk.ac.ox.well.indiana.commands.visualization;

import com.google.common.base.Joiner;
import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultDirectedGraph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.jgrapht.MultiWeightEdge;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

@Description(text="Produce a DAG representation of the contents of a FASTA file")
public class seqgraph extends Module {
    @Argument(fullName="sequences", shortName="s", doc="Sequences FASTA")
    public FastaSequenceFile FASTA;

    @Argument(fullName="highlight", shortName="hl", doc="Colon-separated regex:color list (e.g. 'Pf3D7:#ff0000') specifying the sequence(s) to highlight with the specified color.", required=false)
    public ArrayList<String> HIGHLIGHT;

    @Argument(fullName="contigs", shortName="c", doc="Contigs to add to FASTA graph", required=false)
    public ArrayList<File> CONTIGS;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Argument(fullName="skipSimplification", shortName="ss", doc="Don't simplify the graph")
    public Boolean SKIP_SIMPLIFICATION = false;

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

    @Argument(fullName="numTopEdges", shortName="nte", doc="Number of top (heaviest and longest) edges to list")
    public Integer NUM_TOP_EDGES = 10;

    @Argument(fullName="numVerticesContext", shortName="nvc", doc="Number of vertices to show for context")
    public Integer NUM_VERTICES_CONTEXT = 10;

    @Output(fullName="out", shortName="o", doc="Output file for sequence graph (tip: use a .dot extension)")
    public File out;

    private int numSequences = 0;

    private String joinAttributes(Map<String, Object> m) {
        List<String> attrs = new ArrayList<String>();

        for (String key : m.keySet()) {
            attrs.add(key + "=\"" + m.get(key) + "\"");
        }

        return Joiner.on(" ").join(attrs);
    }

    public void writeGraph(DirectedGraph<String, MultiWeightEdge> g, File f) {
        try {
            PrintStream ps = new PrintStream(f);

            String indent = "  ";

            ps.println("digraph G {");
            ps.println(indent + "rankdir=\"LR\";");

            for (String vertex : g.vertexSet()) {
                boolean isInvisible = true;
                for (MultiWeightEdge e : g.incomingEdgesOf(vertex)) {
                    Map<String, Integer> weights = e.getWeights();

                    for (String color : weights.keySet()) {
                        if (!color.contains("ffffff")) {
                            isInvisible = false;
                            break;
                        }
                    }
                }

                for (MultiWeightEdge e : g.outgoingEdgesOf(vertex)) {
                    Map<String, Integer> weights = e.getWeights();

                    for (String color : weights.keySet()) {
                        if (!color.contains("ffffff")) {
                            isInvisible = false;
                            break;
                        }
                    }
                }

                Map<String, Object> vertexAttrs = new TreeMap<String, Object>();
                vertexAttrs.put("label", WITH_VERTEX_LABELS ? vertex : "");
                vertexAttrs.put("color", "black");
                vertexAttrs.put("shape", WITH_VERTEX_LABELS ? "rect" : "circle");
                vertexAttrs.put("height", 0.10);
                vertexAttrs.put("fontsize", WITH_VERTEX_LABELS ? 10 : 1);
                if (isInvisible) { vertexAttrs.put("style", "invis"); }

                String attributeStr = joinAttributes(vertexAttrs);

                ps.println(indent + vertex + " [ " + attributeStr + " ];");
            }

            for (MultiWeightEdge e : g.edgeSet()) {
                Map<String, Integer> weights = e.getWeights();

                if (weights.size() > 0) {
                    for (String color : weights.keySet()) {
                        float penwidth = MIN_THICKNESS + weights.get(color)*((MAX_THICKNESS - MIN_THICKNESS) / (numSequences - 1));
                        if (Float.isInfinite(penwidth)) {
                            penwidth = MAX_THICKNESS;
                        }

                        Map<String, Object> edgeAttrs = new TreeMap<String, Object>();
                        edgeAttrs.put("penwidth", penwidth);
                        edgeAttrs.put("arrowsize", 0.8);
                        edgeAttrs.put("arrowhead", "normal");
                        edgeAttrs.put("color", color);
                        edgeAttrs.put("samples", Joiner.on(",").join(e.getSequences()));
                        if (color.contains("ffffff")) { edgeAttrs.put("style", "invis"); }
                        if (WITH_EDGE_LABELS) { edgeAttrs.put("label", weights.get(color)); }
                        if (WITH_WEIGHT) { edgeAttrs.put("weight", weights.get(color)); }

                        String attributeStr = joinAttributes(edgeAttrs);

                        ps.println(indent + g.getEdgeSource(e) + " -> " + g.getEdgeTarget(e) + " [ " + attributeStr + " ];");
                    }
                } else {
                    ps.println(indent + g.getEdgeSource(e) + " -> " + g.getEdgeTarget(e));
                }
            }

            ps.println("}");
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }
    }

    private class EdgeComparator implements Comparator<MultiWeightEdge> {
        private DirectedGraph<String, MultiWeightEdge> graph;

        public EdgeComparator(DirectedGraph<String, MultiWeightEdge> graph) {
            this.graph = graph;
        }

        @Override
        public int compare(MultiWeightEdge o1, MultiWeightEdge o2) {
            if (o1.totalWeight() == o2.totalWeight()) {
                int l1 = graph.getEdgeSource(o1).length() + graph.getEdgeTarget(o1).length();
                int l2 = graph.getEdgeSource(o2).length() + graph.getEdgeTarget(o2).length();

                if (l1 == l2) { return 0; }

                return l1 < l2 ? 1 : -1;
            }

            return o1.totalWeight() < o2.totalWeight() ? 1 : -1;
        }
    }

    @Override
    public void execute() {
        DirectedGraph<String, MultiWeightEdge> graph = new DefaultDirectedGraph<String, MultiWeightEdge>(MultiWeightEdge.class);

        // Load genes
        ReferenceSequence rseq;
        Map<String, String> sequences = new HashMap<String, String>();
        while ((rseq = FASTA.nextSequence()) != null) {
            // Load each gene into its own graph object.
            DirectedGraph<String, MultiWeightEdge> sampleGraph = new DefaultDirectedGraph<String, MultiWeightEdge>(MultiWeightEdge.class);
            Set<String> prevKmers = new HashSet<String>();
            String prevKmer = null;

            String seq = new String(rseq.getBases());

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                String kmer = seq.substring(i, i + KMER_SIZE);

                sampleGraph.addVertex(kmer);

                if (prevKmer != null && !prevKmers.contains(kmer)) {
                    sampleGraph.addEdge(prevKmer, kmer);
                    sampleGraph.getEdge(prevKmer, kmer).incrementWeight("#000000");
                }

                prevKmer = kmer;

                if (!SKIP_SIMPLIFICATION) {
                    prevKmers.add(kmer);
                }
            }

            // Add all the graphs together.
            Graphs.addGraph(graph, sampleGraph);

            String[] name = rseq.getName().split("\\s+");

            sequences.put(name[0], seq);
            numSequences++;
        }

        // Load contigs
        if (CONTIGS != null && CONTIGS.size() > 0) {
            DirectedGraph<String, MultiWeightEdge> contigGraph = new DefaultDirectedGraph<String, MultiWeightEdge>(MultiWeightEdge.class);
            for (File contigsFile : CONTIGS) {
                TableReader tr = new TableReader(contigsFile);

                for (Map<String, String> te : tr) {
                    String contigFw = te.get("contig");
                    String contigRc = SequenceUtils.reverseComplement(contigFw);

                    int fwKmersFound = 0;
                    int rcKmersFound = 0;

                    for (int i = 0; i <= contigFw.length() - KMER_SIZE; i++) {
                        String kmerFw = contigFw.substring(i, i + KMER_SIZE);
                        String kmerRc = contigRc.substring(i, i + KMER_SIZE);

                        if (graph.containsVertex(kmerFw)) { fwKmersFound++; }
                        if (graph.containsVertex(kmerRc)) { rcKmersFound++; }
                    }

                    String contig = (fwKmersFound > rcKmersFound) ? contigFw : contigRc;

                    Set<String> prevKmers = new HashSet<String>();
                    String prevKmer = null;

                    for (int i = 0; i <= contig.length() - KMER_SIZE; i++) {
                        String kmer = contig.substring(i, i + KMER_SIZE);

                        contigGraph.addVertex(kmer);

                        if (prevKmer != null && !prevKmers.contains(kmer)) {
                            contigGraph.addEdge(prevKmer, kmer);
                            contigGraph.getEdge(prevKmer, kmer).incrementWeight("#000000");
                        }

                        prevKmer = kmer;

                        if (!SKIP_SIMPLIFICATION) {
                            prevKmers.add(kmer);
                        }
                    }

                    sequences.put(te.get("sample") + ":" + contig.hashCode(), contig);
                }
            }

            Graphs.addGraph(graph, contigGraph);
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

        if (SKIP_SIMPLIFICATION) {
            branchVertices.addAll(graph.vertexSet());
        } else {
            for (String vertex : graph.vertexSet()) {
                if (graph.inDegreeOf(vertex) != 1 || graph.outDegreeOf(vertex) != 1) {
                    branchVertices.add(vertex);
                }
            }
        }

        // Build the complete graph again, but this time, break the sequences into
        // variable sized chunks based on where the branching kmers are found and
        // annotate the edges with weight information.
        DirectedGraph<String, MultiWeightEdge> finalGraph;
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
                        finalGraph.getEdge(prevContig, contig.toString()).addSequence(seqName);
                    }

                    prevContig = contig.toString();
                    contig = new StringBuilder();
                    contig.append(kmer);
                }
            }
        }

        log.info("Sequences: {}", sequences.size());
        log.info("Edges: {}", graph.edgeSet().size());
        log.info("Vertices: {}", graph.vertexSet().size());
        log.info("Branch vertices: {}", branchVertices.size());
        log.info("Simplified vertices: {}", finalGraph.vertexSet().size());
        log.info("Simplified edges: {}", finalGraph.edgeSet().size());

        List<MultiWeightEdge> sortedEdges = new ArrayList<MultiWeightEdge>();
        sortedEdges.addAll(finalGraph.edgeSet());
        Collections.sort(sortedEdges, new EdgeComparator(finalGraph));

        log.info("Top {} heaviest edges with {} vertices of context: ", NUM_TOP_EDGES, NUM_VERTICES_CONTEXT);
        for (int i = 0; i < NUM_TOP_EDGES; i++) {
            log.info("\t{}: {} {}", i, sortedEdges.get(i).totalWeight(), sortedEdges.get(i));
        }

        // Write the graph to disk.
        writeGraph(finalGraph, out);
    }
}
