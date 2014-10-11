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
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexJunctionsRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;
import java.util.List;

public class VisualizeLocalGraph extends Module {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph (sorted)")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="cortexPaths", shortName="cp", doc="Cortex paths")
    public CortexLinks CORTEX_PATHS;

    @Argument(fullName="contigs", shortName="c", doc="Contigs FASTA")
    public IndexedFastaSequenceFile CONTIGS;

    @Argument(fullName="contigName", shortName="cn", doc="Contig name")
    public String CONTIG_NAME;

    @Argument(fullName="showVertexLabels", shortName="vl", doc="Show vertex labels")
    public Boolean SHOW_VERTEX_LABELS = false;

    @Output
    public File out;

    private class WeightedEdge extends DefaultWeightedEdge {
        private String edge;
        private double weight = 1.0;
        private Set<CortexJunctionsRecord> pathWeight = new LinkedHashSet<CortexJunctionsRecord>();

        public WeightedEdge(String edge, double weight) {
            this.edge = edge;
            this.weight = weight;
        }

        public String getLabel() { return edge; }

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

    private void writeGraph(DirectedGraph<String, WeightedEdge> dbg, Set<DirectedGraph<String, WeightedEdge>> pgs, Map<CortexKmer, CortexLinksRecord> paths, Set<String> contigKmers, PrintStream o) {
        String indent = "  ";

        o.println("digraph G {");
        o.println(indent + "rankdir=\"LR\";");

        for (String vertex : dbg.vertexSet()) {
            Map<String, Object> vertexAttrs = new TreeMap<String, Object>();

            CortexKmer ck = new CortexKmer(vertex);
            CortexRecord cr = CORTEX_GRAPH.findRecord(ck);

            if (cr == null) {
                throw new IndianaException("Somehow got to a kmer that doesn't exist in the Cortex binary.");
            }

            String edges = ck.isFlipped() ? SequenceUtils.reverseComplement(cr.getEdgesAsString(0)) : cr.getEdgesAsString(0);
            String label = vertex + "\n" + cr.getCoverage(0) + " " + edges;

            vertexAttrs.put("color", "black");
            vertexAttrs.put("height", 0.10);
            vertexAttrs.put("fontname", "Courier New");
            vertexAttrs.put("fontsize", 8.0);
            if (SHOW_VERTEX_LABELS) {
                vertexAttrs.put("shape", "rect");
                vertexAttrs.put("label", label);
            } else {
                vertexAttrs.put("shape", "circle");
                vertexAttrs.put("label", "");
            }

            if (dbg.outDegreeOf(vertex) > 1 || dbg.inDegreeOf(vertex) > 1 || cr.getInEdgesAsStrings(0).size() > 1 || cr.getOutEdgesAsStrings(0).size() > 1) {
                vertexAttrs.put("shape", "octagon");
                vertexAttrs.put("label", label);
            }

            if (paths.containsKey(ck)) {
                vertexAttrs.put("color", "red");
                if (vertexAttrs.get("shape").equals("circle")) {
                    vertexAttrs.put("shape", "rect");
                }
                vertexAttrs.put("label", label + "\n\n" + (ck.isFlipped() ? reverseComplementJunctions(paths.get(ck)) : paths.get(ck)));
            }

            if (!contigKmers.contains(vertex)) {
                vertexAttrs.put("style", "dashed");
            }

            String attributeStr = joinAttributes(vertexAttrs);

            o.println(indent + vertex + " [ " + attributeStr + " ];");
        }

        for (WeightedEdge e : dbg.edgeSet()) {
            String[] labels = e.getLabel().split("\\n");

            Map<String, Object> edgeAttrs = new TreeMap<String, Object>();
            edgeAttrs.put("arrowsize", 0.8);
            edgeAttrs.put("arrowhead", "normal");
            edgeAttrs.put("penwidth", e.getWeight());
            edgeAttrs.put("weight", e.getWeight());
            edgeAttrs.put("fontname", "Courier New");
            edgeAttrs.put("fontsize", 8.0);
            edgeAttrs.put("dir", "none");
            edgeAttrs.put("color", "black");

            if (!contigKmers.contains(dbg.getEdgeTarget(e)) || !contigKmers.contains(dbg.getEdgeSource(e))) {
                edgeAttrs.put("style", "dashed");
            }

            edgeAttrs.put("label", labels[0]);
            String attributeStr = joinAttributes(edgeAttrs);
            o.println(indent + dbg.getEdgeSource(e) + " -> " + dbg.getEdgeTarget(e) + " [ " + attributeStr + " ];");

            edgeAttrs.put("color", "white");
            edgeAttrs.put("label", labels[1]);
            attributeStr = joinAttributes(edgeAttrs);
            o.println(indent + dbg.getEdgeSource(e) + " -> " + dbg.getEdgeTarget(e) + " [ " + attributeStr + " ];");
        }

        Random rdm = new Random();

        for (DirectedGraph<String, WeightedEdge> pg : pgs) {
            String style = styles[rdm.nextInt(styles.length)];
            String color = colors[rdm.nextInt(colors.length)];

            for (WeightedEdge e : pg.edgeSet()) {
                String label = e.getLabel();

                Map<String, Object> edgeAttrs = new TreeMap<String, Object>();
                edgeAttrs.put("arrowsize", 0.4);
                edgeAttrs.put("arrowhead", "normal");
                edgeAttrs.put("penwidth", e.getWeight());
                edgeAttrs.put("weight", e.getWeight());
                edgeAttrs.put("fontname", "Courier New");
                edgeAttrs.put("fontsize", 8.0);
                edgeAttrs.put("color", color);
                edgeAttrs.put("style", style);
                edgeAttrs.put("label", label);

                String attributeStr = joinAttributes(edgeAttrs);
                o.println(indent + dbg.getEdgeSource(e) + " -> " + dbg.getEdgeTarget(e) + " [ " + attributeStr + " ];");
            }
        }

        o.println("}");
    }

    private Map<CortexKmer, CortexLinksRecord> loadPaths() {
        Map<CortexKmer, CortexLinksRecord> paths = new HashMap<CortexKmer, CortexLinksRecord>();
        for (CortexLinksRecord cpr : CORTEX_PATHS) {
            paths.put(cpr.getKmer(), cpr);
        }

        return paths;
    }

    private String getNextKmer(CortexGraph cg, String kmer) {
        CortexKmer ck = new CortexKmer(kmer);
        CortexRecord cr = cg.findRecord(ck);

        if (cr != null) {
            Collection<String> outEdges = cr.getOutEdgesAsStrings(0);

            if (ck.isFlipped()) {
                outEdges = cr.getInEdgesComplementAsStrings(0);
            }

            if (outEdges.size() == 1) {
                String outEdge = outEdges.iterator().next();

                return kmer.substring(1, kmer.length()) + outEdge;
            }
        }

        return null;
    }

    private Set<String> getNextKmers(CortexGraph cg, String kmer) {
        CortexKmer ck = new CortexKmer(kmer);
        CortexRecord cr = cg.findRecord(ck);
        Set<String> nextKmers = new HashSet<String>();

        if (cr != null) {
            Collection<String> outEdges = cr.getOutEdgesAsStrings(0);

            if (ck.isFlipped()) {
                outEdges = cr.getInEdgesComplementAsStrings(0);
            }


            for (String outEdge : outEdges) {
                nextKmers.add(kmer.substring(1, kmer.length()) + outEdge);
            }
        }

        return nextKmers;
    }

    private String getPrevKmer(CortexGraph cg, String kmer) {
        CortexKmer ck = new CortexKmer(kmer);
        CortexRecord cr = cg.findRecord(ck);

        if (cr != null) {
            Collection<String> inEdges = cr.getInEdgesAsStrings(0);

            if (ck.isFlipped()) {
                inEdges = cr.getOutEdgesComplementAsStrings(0);
            }

            if (inEdges.size() == 1) {
                String inEdge = inEdges.iterator().next();

                return inEdge + kmer.substring(0, kmer.length() - 1);
            }
        }

        return null;
    }

    private Set<String> getPrevKmers(CortexGraph cg, String kmer) {
        CortexKmer ck = new CortexKmer(kmer);
        CortexRecord cr = cg.findRecord(ck);
        Set<String> prevKmers = new HashSet<String>();

        if (cr != null) {
            Collection<String> inEdges = cr.getInEdgesAsStrings(0);

            if (ck.isFlipped()) {
                inEdges = cr.getOutEdgesComplementAsStrings(0);
            }

            for (String inEdge : inEdges) {
                prevKmers.add(inEdge + kmer.substring(0, kmer.length() - 1));
            }
        }

        return prevKmers;
    }

    @Override
    public void execute() {
        int kmerSize = CORTEX_GRAPH.getKmerSize();

        Map<CortexKmer, CortexLinksRecord> paths = loadPaths();

        /*
        if (SHOW_CANDIDATES) {
            ReferenceSequence rseq;
            while ((rseq = CONTIGS.nextSequence()) != null) {
                String seq = new String(rseq.getBases());

                Set<CortexKmer> kmersWithPaths = new HashSet<CortexKmer>();

                for (int i = 0; i <= seq.length() - kmerSize; i++) {
                    String kmer = seq.substring(i, i + kmerSize);

                    CortexKmer ck = new CortexKmer(kmer);
                    if (paths.containsKey(ck)) {
                        kmersWithPaths.add(ck);
                    }
                }

                if (kmersWithPaths.size() == 2 && seq.length() > 100 && seq.length() < 400) {
                    log.info("{}: {} {}", rseq.getName(), seq.length(), kmersWithPaths.size());
                    log.info("  {}", paths.get(kmersWithPaths.iterator().next()));
                }
            }
        }
        */

        log.info("Constructing graph from contigs...");
        String seq = new String(CONTIGS.getSequence(CONTIG_NAME).getBases());
        DirectedGraph<String, WeightedEdge> dbg = new DefaultDirectedGraph<String, WeightedEdge>(WeightedEdge.class);
        String curKmer = null;
        Set<String> contigKmers = new HashSet<String>();

        for (int i = 0; i <= seq.length() - kmerSize; i++) {
            String kmer = seq.substring(i, i + kmerSize);

            contigKmers.add(kmer);
            dbg.addVertex(kmer);

            if (curKmer != null) {
                if (dbg.containsEdge(curKmer, kmer)) {
                    dbg.getEdge(curKmer, kmer).setWeight(dbg.getEdge(curKmer, kmer).getWeight() + 1);
                } else {
                    dbg.addEdge(curKmer, kmer, new WeightedEdge(kmer.substring(kmer.length() - 1, kmer.length()) + "\n" + curKmer.substring(0, 1), 4.0));
                }
            }

            curKmer = kmer;
        }

        log.info("Adding nearby vertices from Cortex graph...");
        int oldVertices = 0;

        for (int i = 0; i < 2; i++) {
            oldVertices = dbg.vertexSet().size();

            log.info("  {} vertices", oldVertices);

            Set<String> vertices = new HashSet<String>();
            vertices.addAll(dbg.vertexSet());

            for (String kmer : vertices) {
                CortexKmer ck = new CortexKmer(kmer);
                CortexRecord cr = CORTEX_GRAPH.findRecord(ck);

                if (cr != null) {
                    Set<String> nextKmers = getNextKmers(CORTEX_GRAPH, kmer);

                    for (String nextKmer : nextKmers) {
                        curKmer = kmer;

                        while (nextKmer != null && !dbg.containsVertex(nextKmer)) {
                            dbg.addVertex(nextKmer);
                            dbg.addEdge(curKmer, nextKmer, new WeightedEdge(nextKmer.substring(nextKmer.length() - 1, nextKmer.length()) + "\n" + curKmer.substring(0, 1), 1.0));

                            curKmer = nextKmer;
                            nextKmer = getNextKmer(CORTEX_GRAPH, nextKmer);
                        }
                    }
                }
            }

            vertices = new HashSet<String>();
            vertices.addAll(dbg.vertexSet());

            for (String kmer : vertices) {
                CortexKmer ck = new CortexKmer(kmer);
                CortexRecord cr = CORTEX_GRAPH.findRecord(ck);

                if (cr != null) {
                    Set<String> prevKmers = getPrevKmers(CORTEX_GRAPH, kmer);

                    for (String prevKmer : prevKmers) {
                        curKmer = kmer;

                        while (prevKmer != null && !dbg.containsVertex(prevKmer)) {
                            dbg.addVertex(prevKmer);
                            dbg.addEdge(prevKmer, curKmer, new WeightedEdge(curKmer.substring(curKmer.length() - 1, curKmer.length()) + "\n" + prevKmer.substring(0, 1), 1.0));

                            curKmer = prevKmer;
                            prevKmer = getPrevKmer(CORTEX_GRAPH, prevKmer);
                        }
                    }
                }
            }
        }

        Set<String> vertices = new HashSet<String>();
        vertices.addAll(dbg.vertexSet());

        Set<DirectedGraph<String, WeightedEdge>> pgs = new HashSet<DirectedGraph<String, WeightedEdge>>();

        log.info("Adding path annotations...");
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
                        while (dbg.containsVertex(curVertex) && junctionsUsed < junctions.length() && dbg.outDegreeOf(curVertex) > 0) {
                            Set<WeightedEdge> nextEdges = dbg.outgoingEdgesOf(curVertex);

                            if (nextEdges.size() == 1) {
                                WeightedEdge edge = nextEdges.iterator().next();
                                String nextVertex = dbg.getEdgeTarget(edge);

                                pg.addVertex(curVertex);
                                pg.addVertex(nextVertex);
                                pg.addEdge(curVertex, nextVertex, new WeightedEdge(junctions.toLowerCase(), 1.0));

                                curVertex = nextVertex;
                            } else {
                                String nextVertex = null;
                                char tBase = 'N';
                                for (WeightedEdge tEdge : nextEdges) {
                                    String tVertex = dbg.getEdgeTarget(tEdge);

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
                                    pg.addEdge(curVertex, nextVertex, new WeightedEdge(junctions.toLowerCase(), 1.0));

                                    curVertex = nextVertex;
                                } else {
                                    log.info("      incomplete traversal at {}: specified junctions did not match available junctions", curVertex);
                                    break;
                                }
                            }
                        }

                        pgs.add(pg);
                    } else {
                        while (dbg.containsVertex(curVertex) && junctionsUsed < junctions.length() && dbg.inDegreeOf(curVertex) > 0) {
                            Set<WeightedEdge> prevEdges = dbg.incomingEdgesOf(curVertex);

                            if (prevEdges.size() == 1) {
                                WeightedEdge edge = prevEdges.iterator().next();
                                String prevVertex = dbg.getEdgeSource(edge);

                                pg.addVertex(curVertex);
                                pg.addVertex(prevVertex);
                                pg.addEdge(curVertex, prevVertex, new WeightedEdge(junctions.toLowerCase(), 1.0));

                                curVertex = prevVertex;
                            } else {
                                String prevVertex = null;
                                char tBase = 'N';
                                for (WeightedEdge tEdge : prevEdges) {
                                    String tVertex = dbg.getEdgeSource(tEdge);

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
                                    pg.addEdge(curVertex, prevVertex, new WeightedEdge(junctions.toLowerCase(), 1.0));

                                    curVertex = prevVertex;
                                } else {
                                    log.info("      incomplete traversal at {}: specified junctions did not match available junctions", curVertex);
                                    break;
                                }
                            }
                        }

                        pgs.add(pg);
                    }
                }
            }
        }

        log.info("Writing graph...");
        try {
            PrintStream o = new PrintStream(out);
            writeGraph(dbg, pgs, paths, contigKmers, o);
        } catch (FileNotFoundException e) {
            throw new IndianaException("Unable to open print stream", e);
        }
    }
}
