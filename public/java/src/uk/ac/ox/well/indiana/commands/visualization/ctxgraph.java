package uk.ac.ox.well.indiana.commands.visualization;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexJunctionsRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksRecord;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

@Description(text="Produce a DAG representation of a region of a Cortex binary file")
public class ctxgraph extends Module {
    @Argument(fullName="inferredEdgesGraph", shortName="ieg", doc="Graph with inferred edges")
    public CortexGraph INFERRED_EDGES_GRAPH;

    @Argument(fullName="uncleanedGraph", shortName="ug", doc="Graph with uncleaned edges", required=false)
    public CortexGraph UNCLEANED_GRAPH;

    @Argument(fullName="paths", shortName="p", doc="McCortex path information", required=false)
    public CortexLinks CORTEX_PATHS;

    @Argument(fullName="contigs", shortName="c", doc="Contigs (in FASTA format)")
    public IndexedFastaSequenceFile CONTIGS;

    @Argument(fullName="contigName", shortName="cn", doc="Contig name", required=false)
    public String CONTIG_NAME;

    @Output
    public PrintStream out;

    private class WeightedEdge extends DefaultWeightedEdge {
        private String inEdge;
        private String outEdge;
        private double weight = 1.0;
        private Set<CortexJunctionsRecord> pathWeight = new LinkedHashSet<CortexJunctionsRecord>();

        public WeightedEdge(String inEdge, String outEdge, double weight) {
            this.inEdge = inEdge;
            this.outEdge = outEdge;
            this.weight = weight;
        }

        public String getInEdge() { return inEdge; }
        public String getOutEdge() { return outEdge; }

        public void setWeight(double weight) { this.weight = weight; }
        public double getWeight() { return weight; }

        public void addPathWeight(CortexJunctionsRecord jr) { pathWeight.add(jr); }
        public int numPaths() { return pathWeight.size(); }
        public Set<CortexJunctionsRecord> getPathWeights() { return pathWeight; }
    }

    private String loadContig() {
        ReferenceSequence rseq = (CONTIG_NAME == null) ? CONTIGS.nextSequence() : CONTIGS.getSequence(CONTIG_NAME);

        log.info("  contig: {}", rseq.getName());

        return new String(rseq.getBases());
    }

    private DirectedGraph<String, WeightedEdge> loadLocalGraph(Set<String> contigKmers) {
        DirectedGraph<String, WeightedEdge> g = new DefaultDirectedGraph<String, WeightedEdge>(WeightedEdge.class);

        for (String contigKmer : contigKmers) {
            CortexKmer ck = new CortexKmer(contigKmer);

            CortexRecord cr = INFERRED_EDGES_GRAPH.findRecord(ck);

            if (cr != null) {
                g.addVertex(contigKmer);

                Set<String> nextKmers = CortexUtils.getNextKmers(INFERRED_EDGES_GRAPH, contigKmer, 0);
                for (String nextKmer : nextKmers) {
                    g.addVertex(nextKmer);

                    String inEdge = contigKmer.substring(0, 1);
                    String outEdge = nextKmer.substring(nextKmer.length() - 1, nextKmer.length());

                    g.addEdge(contigKmer, nextKmer, new WeightedEdge(inEdge, outEdge, 1.0));
                }

                Set<String> prevKmers = CortexUtils.getPrevKmers(INFERRED_EDGES_GRAPH, contigKmer, 0);
                for (String prevKmer : prevKmers) {
                    g.addVertex(prevKmer);

                    String inEdge = prevKmer.substring(0, 1);
                    String outEdge = contigKmer.substring(contigKmer.length() - 1, contigKmer.length());

                    g.addEdge(prevKmer, contigKmer, new WeightedEdge(inEdge, outEdge, 1.0));
                }
            }
        }

        return g;
    }

    private void extendLocalGraph(DirectedGraph<String, WeightedEdge> g, int iterations) {
        int oldVertices = 0;

        for (int i = 0; i < iterations; i++) {
            oldVertices = g.vertexSet().size();

            log.info("  {} vertices", oldVertices);

            Set<String> vertices = new HashSet<String>();
            vertices.addAll(g.vertexSet());

            String curKmer = null;
            for (String kmer : vertices) {
                CortexKmer ck = new CortexKmer(kmer);
                CortexRecord cr = INFERRED_EDGES_GRAPH.findRecord(ck);

                if (cr != null) {
                    Set<String> nextKmers = CortexUtils.getNextKmers(INFERRED_EDGES_GRAPH, kmer, 0);

                    for (String nextKmer : nextKmers) {
                        curKmer = kmer;

                        while (nextKmer != null && !g.containsVertex(nextKmer)) {
                            g.addVertex(nextKmer);

                            String outEdge = nextKmer.substring(nextKmer.length() - 1, nextKmer.length());
                            String inEdge = curKmer.substring(0, 1);

                            g.addEdge(curKmer, nextKmer, new WeightedEdge(inEdge, outEdge, 1.0));

                            curKmer = nextKmer;
                            nextKmer = CortexUtils.getNextKmer(INFERRED_EDGES_GRAPH, nextKmer, 0, false);
                        }
                    }
                }
            }

            vertices = new HashSet<String>();
            vertices.addAll(g.vertexSet());

            for (String kmer : vertices) {
                CortexKmer ck = new CortexKmer(kmer);
                CortexRecord cr = INFERRED_EDGES_GRAPH.findRecord(ck);

                if (cr != null) {
                    Set<String> prevKmers = CortexUtils.getPrevKmers(INFERRED_EDGES_GRAPH, kmer, 0);

                    for (String prevKmer : prevKmers) {
                        curKmer = kmer;

                        while (prevKmer != null && !g.containsVertex(prevKmer)) {
                            g.addVertex(prevKmer);

                            String inEdge = prevKmer.substring(0, 1);
                            String outEdge = curKmer.substring(curKmer.length() - 1, curKmer.length());

                            g.addEdge(prevKmer, curKmer, new WeightedEdge(inEdge, outEdge, 1.0));

                            curKmer = prevKmer;
                            prevKmer = CortexUtils.getPrevKmer(INFERRED_EDGES_GRAPH, prevKmer, 0, false);
                        }
                    }
                }
            }
        }
    }

    private Map<CortexKmer, CortexLinksRecord> loadPaths() {
        Map<CortexKmer, CortexLinksRecord> paths = new HashMap<CortexKmer, CortexLinksRecord>();

        if (CORTEX_PATHS != null) {
            for (CortexLinksRecord cpr : CORTEX_PATHS) {
                paths.put(cpr.getKmer(), cpr);
            }
        }

        return paths;
    }

    private Set<DirectedGraph<String, WeightedEdge>> addPathAnnotations(DirectedGraph<String, WeightedEdge> g, Map<CortexKmer, CortexLinksRecord> paths) {
        Set<String> vertices = new HashSet<String>();
        vertices.addAll(g.vertexSet());

        Set<DirectedGraph<String, WeightedEdge>> pgs = new HashSet<DirectedGraph<String, WeightedEdge>>();

        for (String vertex : vertices) {
            CortexKmer ck = new CortexKmer(vertex);

            if (paths.containsKey(ck)) {
                CortexLinksRecord pr = paths.get(ck);
                log.info("  paths starting at: {}", vertex);

                for (CortexJunctionsRecord jr : pr.getJunctions()) {
                    log.info("    {}", jr);

                    DirectedGraph<String, WeightedEdge> pg = new DefaultDirectedGraph<String, WeightedEdge>(WeightedEdge.class);

                    boolean goForward = jr.isForward();
                    String junctions = jr.getJunctions();
                    int junctionsUsed = 0;

                    if (ck.isFlipped()) {
                        goForward = !goForward;
                        junctions = SequenceUtils.complement(junctions);
                    }

                    String curVertex = vertex;
                    if (goForward) {
                        boolean completeTraversal = true;

                        while (g.containsVertex(curVertex) && junctionsUsed < junctions.length() && g.outDegreeOf(curVertex) > 0) {
                            Set<WeightedEdge> nextEdges = g.outgoingEdgesOf(curVertex);

                            if (nextEdges.size() == 1) {
                                WeightedEdge edge = nextEdges.iterator().next();
                                String nextVertex = g.getEdgeTarget(edge);

                                pg.addVertex(curVertex);
                                pg.addVertex(nextVertex);
                                pg.addEdge(curVertex, nextVertex, new WeightedEdge("", junctions.toLowerCase(), 1.0));

                                curVertex = nextVertex;
                            } else {
                                String nextVertex = null;
                                char tBase = 'N';
                                for (WeightedEdge tEdge : nextEdges) {
                                    String tVertex = g.getEdgeTarget(tEdge);

                                    tBase = tVertex.substring(tVertex.length() - 1, tVertex.length()).charAt(0);

                                    if (tBase == junctions.charAt(junctionsUsed)) {
                                        nextVertex = tVertex;
                                        junctionsUsed++;
                                        break;
                                    }
                                }

                                if (nextVertex != null) {
                                    pg.addVertex(curVertex);
                                    pg.addVertex(nextVertex);
                                    pg.addEdge(curVertex, nextVertex, new WeightedEdge("", junctions.toLowerCase(), 1.0));

                                    curVertex = nextVertex;
                                } else {
                                    log.info("      incomplete traversal at {}: specified junctions did not match available junctions", curVertex);
                                    completeTraversal = false;
                                    break;
                                }
                            }
                        }

                        if (completeTraversal) {
                            pgs.add(pg);
                        }
                    } else {
                        boolean completeTraversal = true;

                        while (g.containsVertex(curVertex) && junctionsUsed < junctions.length() && g.inDegreeOf(curVertex) > 0) {
                            Set<WeightedEdge> prevEdges = g.incomingEdgesOf(curVertex);

                            if (prevEdges.size() == 1) {
                                WeightedEdge edge = prevEdges.iterator().next();
                                String prevVertex = g.getEdgeSource(edge);

                                pg.addVertex(curVertex);
                                pg.addVertex(prevVertex);
                                pg.addEdge(curVertex, prevVertex, new WeightedEdge("", junctions.toLowerCase(), 1.0));

                                curVertex = prevVertex;
                            } else {
                                String prevVertex = null;
                                char tBase = 'N';
                                for (WeightedEdge tEdge : prevEdges) {
                                    String tVertex = g.getEdgeSource(tEdge);

                                    tBase = tVertex.substring(0, 1).charAt(0);

                                    if (tBase == junctions.charAt(junctionsUsed)) {
                                        prevVertex = tVertex;
                                        junctionsUsed++;
                                        break;
                                    }
                                }

                                if (prevVertex != null) {
                                    pg.addVertex(curVertex);
                                    pg.addVertex(prevVertex);
                                    pg.addEdge(curVertex, prevVertex, new WeightedEdge("", junctions.toLowerCase(), 1.0));

                                    curVertex = prevVertex;
                                } else {
                                    log.info("      incomplete traversal at {}: specified junctions did not match available junctions", curVertex);
                                    completeTraversal = false;
                                    break;
                                }
                            }
                        }

                        if (completeTraversal) {
                            pgs.add(pg);
                        }
                    }
                }
            }
        }

        return pgs;
    }

    private void addUncleanedKmers(DirectedGraph<String, WeightedEdge> g) {
        if (UNCLEANED_GRAPH != null) {
            Set<String> vertices = new HashSet<String>();
            vertices.addAll(g.vertexSet());

            for (String vertex : vertices) {
                CortexKmer ck = new CortexKmer(vertex);
                CortexRecord cr = UNCLEANED_GRAPH.findRecord(ck);

                if (cr != null) {
                    g.addVertex(vertex);

                    Set<String> nextKmers = CortexUtils.getNextKmers(UNCLEANED_GRAPH, vertex, 0);
                    for (String nextKmer : nextKmers) {
                        g.addVertex(nextKmer);

                        String inEdge = vertex.substring(0, 1);
                        String outEdge = nextKmer.substring(nextKmer.length() - 1, nextKmer.length());

                        g.addEdge(vertex, nextKmer, new WeightedEdge(inEdge, outEdge, 1.0));
                    }

                    Set<String> prevKmers = CortexUtils.getPrevKmers(UNCLEANED_GRAPH, vertex, 0);
                    for (String prevKmer : prevKmers) {
                        g.addVertex(prevKmer);

                        String inEdge = prevKmer.substring(0, 1);
                        String outEdge = vertex.substring(vertex.length() - 1, vertex.length());

                        g.addEdge(prevKmer, vertex, new WeightedEdge(inEdge, outEdge, 1.0));
                    }
                }
            }
        }
    }

    private void addContigKmers(DirectedGraph<String, WeightedEdge> g, String contig, int kmerSize) {
        String prevKmer = null;
        for (int i = 0; i <= contig.length() - kmerSize; i++) {
            String curKmer = contig.substring(i, i + kmerSize);

            g.addVertex(curKmer);

            if (prevKmer != null) {
                String inEdge = prevKmer.substring(0, 1);
                String outEdge = curKmer.substring(curKmer.length() - 1, curKmer.length());

                if (g.containsEdge(prevKmer, curKmer)) {
                    g.getEdge(prevKmer, curKmer).setWeight(5.0);
                } else {
                    g.addEdge(prevKmer, curKmer, new WeightedEdge(inEdge, outEdge, 5.0));
                }
            }

            prevKmer = curKmer;
        }
    }

    private String joinAttributes(Map<String, Object> m) {
        List<String> attrs = new ArrayList<String>();

        for (String key : m.keySet()) {
            attrs.add(key + "=\"" + m.get(key) + "\"");
        }

        return Joiner.on(" ").join(attrs);
    }

    private String reverseComplementJunction(CortexJunctionsRecord jr) {
        StringBuilder buffer = new StringBuilder();

        buffer.append(jr.isForward() ? "R" : "F").append(" ");
        buffer.append(jr.getNumKmers()).append(" ");
        buffer.append(jr.getNumJunctions()).append(" ");

        for (int c = 0; c < jr.getCoverages().length; c++) {
            buffer.append(jr.getCoverage(c)).append(" ");
        }

        buffer.append(SequenceUtils.complement(jr.getJunctions()));

        return buffer.toString();
    }

    private String reverseComplementJunctions(CortexLinksRecord pr) {
        StringBuilder record = new StringBuilder();

        record.append(SequenceUtils.reverseComplement(pr.getKmerAsString())).append(" ").append(pr.getJunctions().size()).append("\n");

        int i = 0;
        for (CortexJunctionsRecord cj : pr.getJunctions()) {
            i++;

            record.append(reverseComplementJunction(cj));

            if (i < pr.getJunctions().size()) { record.append("\n"); }
        }

        return record.toString();
    }

    private void writeGraph(DirectedGraph<String, WeightedEdge> g, Set<DirectedGraph<String, WeightedEdge>> pgs, Map<CortexKmer, CortexLinksRecord> paths, Set<String> contigKmers) {
        Map<String, Object> globalVertexAttrs = new HashMap<String, Object>();
        globalVertexAttrs.put("color", "black");
        globalVertexAttrs.put("fontname", "Courier New");
        globalVertexAttrs.put("fontsize", 7.0);
        globalVertexAttrs.put("height", 0.10);
        globalVertexAttrs.put("shape", "circle");
        globalVertexAttrs.put("label", "");
        globalVertexAttrs.put("style", "dashed");

        Map<String, Object> globalEdgeAttrs = new HashMap<String, Object>();
        globalEdgeAttrs.put("arrowsize", 0.4);
        globalEdgeAttrs.put("arrowhead", "normal");
        globalEdgeAttrs.put("penwidth", 1.0);
        globalEdgeAttrs.put("weight", 1.0);
        globalEdgeAttrs.put("fontname", "Courier New");
        globalEdgeAttrs.put("fontsize", 7.0);
        globalEdgeAttrs.put("color", "black");
        globalEdgeAttrs.put("style", "dashed");

        final String[] colors = { "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6" };
        final String[] styles = { "dashed", "dotted", "solid", "bold" };

        String indent = "  ";

        out.println("digraph G {");
        out.println(indent + "rankdir=\"LR\";");
        out.println(indent + "node [" + joinAttributes(globalVertexAttrs) + "];");
        out.println(indent + "edge [" + joinAttributes(globalEdgeAttrs) + "];");

        for (String vertex : g.vertexSet()) {
            Map<String, Object> vertexAttrs = new HashMap<String, Object>();

            if (contigKmers.contains(vertex)) {
                vertexAttrs.put("style", "solid");
            }

            if (g.inDegreeOf(vertex) == 0 || g.outDegreeOf(vertex) == 0) {
                vertexAttrs.put("label", "...");
                vertexAttrs.put("fontname", "Courier New italic");
            }

            if (g.inDegreeOf(vertex) > 1 || g.outDegreeOf(vertex) > 1) {
                vertexAttrs.put("label", vertex);
                vertexAttrs.put("shape", "octagon");
            }

            CortexKmer ck = new CortexKmer(vertex);
            if (paths.containsKey(ck)) {
                vertexAttrs.put("label", (ck.isFlipped() ? reverseComplementJunctions(paths.get(ck)) : paths.get(ck)));
                vertexAttrs.put("color", "red");
                if (!vertex.contains("shape")) {
                    vertexAttrs.put("shape", "rect");
                }
            }

            if (vertexAttrs.size() > 0) {
                out.println(indent + vertex + " [" + joinAttributes(vertexAttrs) + "];");
            } else {
                out.println(indent + vertex);
            }
        }

        for (WeightedEdge e : g.edgeSet()) {
            String vertexSource = g.getEdgeSource(e);
            String vertexTarget = g.getEdgeTarget(e);

            Map<String, Object> edgeAttrs = new HashMap<String, Object>();
            edgeAttrs.put("label", e.getOutEdge());
            edgeAttrs.put("xlabel", e.getInEdge());

            if (contigKmers.contains(vertexSource) && contigKmers.contains(vertexTarget)) {
                edgeAttrs.put("style", "solid");
                edgeAttrs.put("penwidth", e.getWeight());
                edgeAttrs.put("weight", e.getWeight());
            }

            out.println(indent + vertexSource + " -> " + vertexTarget + " [" + joinAttributes(edgeAttrs) + "];");
        }

        Random rdm = new Random();

        for (DirectedGraph<String, WeightedEdge> pg : pgs) {
            String style = styles[rdm.nextInt(styles.length)];
            String color = colors[rdm.nextInt(colors.length)];

            for (WeightedEdge e : pg.edgeSet()) {
                String label = e.getOutEdge();

                Map<String, Object> edgeAttrs = new TreeMap<String, Object>();
                edgeAttrs.put("penwidth", 0.5);
                edgeAttrs.put("weight", 0.0);
                edgeAttrs.put("color", color);
                edgeAttrs.put("fontcolor", color);
                edgeAttrs.put("style", style);
                edgeAttrs.put("label", label);
                edgeAttrs.put("fontsize", 3.0);

                out.println(indent + g.getEdgeSource(e) + " -> " + g.getEdgeTarget(e) + " [ " + joinAttributes(edgeAttrs) + " ];");
            }
        }

        out.println("}");
    }

    @Override
    public void execute() {
        int kmerSize = INFERRED_EDGES_GRAPH.getKmerSize();

        log.info("Loading selected contig...");
        String contig = loadContig();
        Set<String> contigKmers = SequenceUtils.kmerizeSequence(contig, kmerSize);

        log.info("Loading paths...");
        Map<CortexKmer, CortexLinksRecord> paths = loadPaths();

        log.info("Constructing initial graph...");
        DirectedGraph<String, WeightedEdge> g = loadLocalGraph(contigKmers);

        log.info("Extending graph with nearby kmers...");
        extendLocalGraph(g, 1);

        log.info("Adding path annotations...");
        Set<DirectedGraph<String, WeightedEdge>> pgs = addPathAnnotations(g, paths);

        log.info("Adding kmers from uncleaned graph...");
        addUncleanedKmers(g);

        log.info("Adding contig kmers...");
        addContigKmers(g, contig, kmerSize);

        log.info("Writing graph...");
        writeGraph(g, pgs, paths, contigKmers);
    }
}
