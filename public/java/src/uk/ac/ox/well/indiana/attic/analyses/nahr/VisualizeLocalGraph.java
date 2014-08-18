package uk.ac.ox.well.indiana.attic.analyses.nahr;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.paths.CortexJunctionsRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.paths.CortexPaths;
import uk.ac.ox.well.indiana.utils.io.cortex.paths.CortexPathsRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;
import java.util.List;

public class VisualizeLocalGraph extends Module {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="cortexPaths", shortName="cp", doc="Cortex paths")
    public CortexPaths CORTEX_PATHS;

    @Argument(fullName="contigs", shortName="c", doc="Contigs FASTA")
    public IndexedFastaSequenceFile CONTIGS;

    @Argument(fullName="contigName", shortName="cn", doc="Contig name")
    public String CONTIG_NAME;

    @Output
    public PrintStream out;

    private class WeightedEdge extends DefaultWeightedEdge {
        private double weight = 1.0;
        private Set<CortexJunctionsRecord> pathWeight = new LinkedHashSet<CortexJunctionsRecord>();

        public WeightedEdge(double weight) { this.weight = weight; }
        public void setWeight(double weight) { this.weight = weight; }
        public double getWeight() { return weight; }

        public void addPathWeight(CortexJunctionsRecord jr) { pathWeight.add(jr); }
        public int numPaths() { return pathWeight.size(); }
        public Set<CortexJunctionsRecord> getPathWeights() { return pathWeight; }
    }

    private final String[] colors = { "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6" };
    private final String[] styles = { "dashed", "dotted", "solid", "bold" };

    private String joinAttributes(Map<String, Object> m) {
        List<String> attrs = new ArrayList<String>();

        for (String key : m.keySet()) {
            attrs.add(key + "=\"" + m.get(key) + "\"");
        }

        return Joiner.on(" ").join(attrs);
    }

    private void writeGraph(DirectedGraph<String, WeightedEdge> g, Map<CortexKmer, CortexPathsRecord> paths) {
        String indent = "  ";

        out.println("digraph G {");
        out.println(indent + "rankdir=\"LR\";");

        for (String vertex : g.vertexSet()) {
            Map<String, Object> vertexAttrs = new TreeMap<String, Object>();

            String label = vertex.charAt(0) + ".." + vertex.charAt(vertex.length() - 1);

            vertexAttrs.put("color", "black");
            vertexAttrs.put("height", 0.10);
            vertexAttrs.put("shape", "rect");
            vertexAttrs.put("label", label + "\n" + SequenceUtils.reverseComplement(label));
            vertexAttrs.put("fontname", "Courier New");

            if (g.outDegreeOf(vertex) > 1 || g.inDegreeOf(vertex) > 1) {
                vertexAttrs.put("shape", "hexagon");
                vertexAttrs.put("height", 0.50);
            }

            if (paths.containsKey(new CortexKmer(vertex))) {
                vertexAttrs.put("color", "red");
                //vertexAttrs.put("label", vertex + "\n" + SequenceUtils.reverseComplement(vertex));
                vertexAttrs.put("label", vertex + "\n" + paths.get(new CortexKmer(vertex)).toString());
            }

            String attributeStr = joinAttributes(vertexAttrs);

            out.println(indent + vertex + " [ " + attributeStr + " ];");
        }

        Map<CortexJunctionsRecord, String> colorMap = new HashMap<CortexJunctionsRecord, String>();
        Map<CortexJunctionsRecord, String> styleMap = new HashMap<CortexJunctionsRecord, String>();
        Random rdm = new Random();

        for (WeightedEdge e : g.edgeSet()) {
            Map<String, Object> edgeAttrs = new TreeMap<String, Object>();
            edgeAttrs.put("arrowsize", 0.8);
            edgeAttrs.put("arrowhead", "normal");
            edgeAttrs.put("color", "black");
            edgeAttrs.put("penwidth", e.getWeight());
            edgeAttrs.put("weight", e.getWeight());

            String attributeStr = joinAttributes(edgeAttrs);
            out.println(indent + g.getEdgeSource(e) + " -> " + g.getEdgeTarget(e) + " [ " + attributeStr + " ];");

            if (e.numPaths() > 0) {
                for (CortexJunctionsRecord jr : e.getPathWeights()) {
                    if (!colorMap.containsKey(jr)) {
                        colorMap.put(jr, colors[rdm.nextInt(colors.length)]);
                    }

                    if (!styleMap.containsKey(jr)) {
                        styleMap.put(jr, styles[rdm.nextInt(styles.length)]);
                    }

                    edgeAttrs.put("color", colorMap.get(jr));
                    edgeAttrs.put("style", styleMap.get(jr));
                    edgeAttrs.put("weight", jr.getCoverage(0));
                    edgeAttrs.put("label", jr.getJunctions());

                    attributeStr = joinAttributes(edgeAttrs);

                    if (jr.isForward()) {
                        out.println(indent + g.getEdgeSource(e) + " -> " + g.getEdgeTarget(e) + " [ " + attributeStr + " ];");
                    } else {
                        out.println(indent + g.getEdgeTarget(e) + " -> " + g.getEdgeSource(e) + " [ " + attributeStr + " ];");
                    }
                }
            }
        }

        out.println("}");
    }

    private Map<CortexKmer, CortexPathsRecord> loadPaths() {
        Map<CortexKmer, CortexPathsRecord> paths = new HashMap<CortexKmer, CortexPathsRecord>();
        for (CortexPathsRecord cpr : CORTEX_PATHS) {
            paths.put(cpr.getKmer(), cpr);
        }

        return paths;
    }

    @Override
    public void execute() {
        String seq = new String(CONTIGS.getSequence(CONTIG_NAME).getBases());
        int kmerSize = CORTEX_GRAPH.getKmerSize();

        Map<CortexKmer, CortexPathsRecord> paths = loadPaths();

        DirectedGraph<String, WeightedEdge> dbg = new DefaultDirectedGraph<String, WeightedEdge>(WeightedEdge.class);

        String lastKmer = null;
        for (int i = 0; i <= seq.length() - kmerSize; i++) {
            String kmer = seq.substring(i, i + kmerSize);
            CortexKmer ck = new CortexKmer(kmer);
            CortexRecord cr = CORTEX_GRAPH.findRecord(ck);

            dbg.addVertex(kmer);

            if (cr != null) {
                Collection<String> edgesIn = ck.isFlipped() ? cr.getOutEdgesComplementAsStrings(0) : cr.getInEdgesAsStrings(0);
                Collection<String> edgesOut = ck.isFlipped() ? cr.getInEdgesComplementAsStrings(0) : cr.getOutEdgesAsStrings(0);

                for (String edgeOut : edgesOut) {
                    String nextKmer = kmer.substring(1, kmer.length()) + edgeOut;

                    dbg.addVertex(nextKmer);
                    dbg.addEdge(kmer, nextKmer, new WeightedEdge(1.0));
                }

                for (String edgeIn : edgesIn) {
                    String prevKmer = edgeIn + kmer.substring(0, kmer.length() - 1);

                    dbg.addVertex(prevKmer);
                    dbg.addEdge(prevKmer, kmer, new WeightedEdge(1.0));
                }
            }

            if (lastKmer != null) {
                if (dbg.containsEdge(lastKmer, kmer)) {
                    dbg.getEdge(lastKmer, kmer).setWeight(5.0);
                } else {
                    dbg.addEdge(lastKmer, kmer, new WeightedEdge(5.0));
                }
            }

            lastKmer = kmer;
        }

        for (String vertex : dbg.vertexSet()) {
            CortexKmer ck = new CortexKmer(vertex);

            if (paths.containsKey(ck)) {
                CortexRecord cr = CORTEX_GRAPH.findRecord(ck);

                for (CortexJunctionsRecord jr : paths.get(ck).getJunctions()) {
                    String junctionVertex = null;

                    if (!ck.isFlipped()) {
                        if (jr.isForward()) {
                            Collection<String> edges = cr.getOutEdgesAsStrings(0);
                            if (edges.size() == 1) { junctionVertex = vertex.substring(1, vertex.length()) + edges.iterator().next(); }
                        } else {
                            Collection<String> edges = cr.getInEdgesAsStrings(0);
                            if (edges.size() == 1) { junctionVertex = edges.iterator().next() + vertex.substring(0, vertex.length() - 1); }
                        }
                    } else {
                        if (jr.isForward()) {
                            Collection<String> edges = cr.getOutEdgesComplementAsStrings(0);
                            if (edges.size() == 1) { junctionVertex = edges.iterator().next() + vertex.substring(0, vertex.length() - 1); }
                        } else {
                            Collection<String> edges = cr.getInEdgesComplementAsStrings(0);
                            if (edges.size() == 1) { junctionVertex = vertex.substring(1, vertex.length()) + edges.iterator().next(); }
                        }
                    }

                    if (junctionVertex == null) {
                        throw new IndianaException("Could not find the junction vertex adjacent to '" + vertex + "'");
                    }

                    log.info("{} -> {} ({} {} {})", vertex, junctionVertex, ck.isFlipped(), jr.isForward(), dbg.containsVertex(junctionVertex));

                    if (vertex.equals("CTCTAATTGATCATTCACTAATTGATCATTCACTAATTGATCATTCT")) {
                        log.info("Found problematic vertex");
                    }

                    String junctions = jr.getJunctions();
                    int usedJunctions = 0;
                    String currentVertex = junctionVertex;

                    while (usedJunctions < jr.getNumJunctions() && dbg.containsVertex(currentVertex)) {
                        CortexKmer currentKmer = new CortexKmer(currentVertex);

                        if (jr.isForward()) {
                            Set<WeightedEdge> edges = dbg.outgoingEdgesOf(currentVertex);

                            if (edges.size() > 0) {
                                String nextVertex;

                                if (edges.size() == 1) {
                                    nextVertex = dbg.getEdgeTarget(edges.iterator().next());
                                } else {
                                    nextVertex = currentVertex.substring(1, currentVertex.length()) + junctions.charAt(usedJunctions);
                                    usedJunctions++;
                                }

                                WeightedEdge supportedEdge = dbg.getEdge(currentVertex, nextVertex);
                                if (supportedEdge != null) {
                                    supportedEdge.addPathWeight(jr);
                                }

                                currentVertex = nextVertex;
                            }
                        } else {
                            Set<WeightedEdge> edges = dbg.incomingEdgesOf(currentVertex);

                            if (edges.size() > 0) {
                                String prevVertex;

                                if (edges.size() == 1) {
                                    prevVertex = dbg.getEdgeSource(edges.iterator().next());
                                } else {
                                    if (currentKmer.isFlipped()) {
                                        prevVertex = SequenceUtils.complement(junctions.charAt(usedJunctions)) + currentVertex.substring(0, currentVertex.length() - 1);
                                    } else {
                                        prevVertex = junctions.charAt(usedJunctions) + currentVertex.substring(0, currentVertex.length() - 1);
                                    }
                                    usedJunctions++;
                                }

                                WeightedEdge supportedEdge = dbg.getEdge(prevVertex, currentVertex);
                                if (supportedEdge != null) {
                                    supportedEdge.addPathWeight(jr);
                                }

                                currentVertex = prevVertex;
                            } else {
                                usedJunctions++;
                            }
                        }
                    }
                }
            }
        }

        writeGraph(dbg, paths);
    }
}
