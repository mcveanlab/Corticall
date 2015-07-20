package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.DirectedGraph;
import org.jgrapht.GraphPath;
import org.jgrapht.Graphs;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.alg.DijkstraShortestPath;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.EdgeReversedGraph;
import org.jgrapht.traverse.DepthFirstIterator;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

public class GenotypeGraph extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="graphRaw", shortName="gr", doc="Graph (raw)", required=false)
    public CortexGraph GRAPH_RAW;

    @Argument(fullName="novelGraph", shortName="n", doc="Graph of novel kmers")
    public CortexGraph NOVEL;

    @Argument(fullName="bed", shortName="b", doc="Bed file describing variants")
    public File BED;

    @Argument(fullName="novelKmerMap", shortName="m", doc="Novel kmer map")
    public File NOVEL_KMER_MAP;

    @Output(fullName="gout", shortName="go", doc="Graph out")
    public File gout;

    @Output
    public PrintStream out;

    private class VariantInfo {
        public String variantId;

        public String vchr;
        public int vstart;
        public int vstop;

        public String type;
        public String denovo;
        public String nahr;
        public int gcindex;
        public String ref;
        public String alt;
        public String leftFlank;
        public String rightFlank;

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            VariantInfo that = (VariantInfo) o;

            if (vstart != that.vstart) return false;
            if (vstop != that.vstop) return false;
            if (gcindex != that.gcindex) return false;
            if (variantId != null ? !variantId.equals(that.variantId) : that.variantId != null) return false;
            if (vchr != null ? !vchr.equals(that.vchr) : that.vchr != null) return false;
            if (type != null ? !type.equals(that.type) : that.type != null) return false;
            if (denovo != null ? !denovo.equals(that.denovo) : that.denovo != null) return false;
            if (nahr != null ? !nahr.equals(that.nahr) : that.nahr != null) return false;
            if (ref != null ? !ref.equals(that.ref) : that.ref != null) return false;
            if (alt != null ? !alt.equals(that.alt) : that.alt != null) return false;
            if (leftFlank != null ? !leftFlank.equals(that.leftFlank) : that.leftFlank != null) return false;
            return !(rightFlank != null ? !rightFlank.equals(that.rightFlank) : that.rightFlank != null);

        }

        @Override
        public int hashCode() {
            int result = variantId != null ? variantId.hashCode() : 0;
            result = 31 * result + (vchr != null ? vchr.hashCode() : 0);
            result = 31 * result + vstart;
            result = 31 * result + vstop;
            result = 31 * result + (type != null ? type.hashCode() : 0);
            result = 31 * result + (denovo != null ? denovo.hashCode() : 0);
            result = 31 * result + (nahr != null ? nahr.hashCode() : 0);
            result = 31 * result + gcindex;
            result = 31 * result + (ref != null ? ref.hashCode() : 0);
            result = 31 * result + (alt != null ? alt.hashCode() : 0);
            result = 31 * result + (leftFlank != null ? leftFlank.hashCode() : 0);
            result = 31 * result + (rightFlank != null ? rightFlank.hashCode() : 0);
            return result;
        }

        @Override
        public String toString() {
            return "VariantInfo{" +
                    "variantId='" + variantId + '\'' +
                    ", vchr='" + vchr + '\'' +
                    ", vstart=" + vstart +
                    ", vstop=" + vstop +
                    ", type='" + type + '\'' +
                    ", denovo='" + denovo + '\'' +
                    ", nahr='" + nahr + '\'' +
                    ", gcindex=" + gcindex +
                    ", ref='" + ref + '\'' +
                    ", alt='" + alt + '\'' +
                    ", leftFlank='" + leftFlank + '\'' +
                    ", rightFlank='" + rightFlank + '\'' +
                    '}';
        }
    }

    private Map<CortexKmer, VariantInfo> loadNovelKmerMap() {
        Map<String, VariantInfo> allVariants = new HashMap<String, VariantInfo>();

        TableReader bed = new TableReader(BED, "vchr", "vstart", "vstop", "info");

        for (Map<String, String> be : bed) {
            VariantInfo vi = new VariantInfo();

            vi.vchr = be.get("vchr");
            vi.vstart = Integer.valueOf(be.get("vstart"));
            vi.vstop = Integer.valueOf(be.get("vstop"));

            String[] kvpairs = be.get("info").split(";");
            for (String kvpair : kvpairs) {
                String[] kv = kvpair.split("=");

                if (kv[0].equals("id")) { vi.variantId = kv[1]; }
                if (kv[0].equals("type")) { vi.type = kv[1]; }
                if (kv[0].equals("denovo")) { vi.denovo = kv[1]; }
                if (kv[0].equals("nahr")) { vi.nahr = kv[1]; }
                if (kv[0].equals("gcindex")) { vi.gcindex = Integer.valueOf(kv[1]); }
                if (kv[0].equals("ref")) { vi.ref = kv[1]; }
                if (kv[0].equals("alt")) { vi.alt = kv[1]; }
                if (kv[0].equals("left")) { vi.leftFlank = kv[1]; }
                if (kv[0].equals("right")) { vi.rightFlank = kv[1]; }
            }

            allVariants.put(vi.variantId, vi);
        }

        Map<CortexKmer, VariantInfo> vis = new HashMap<CortexKmer, VariantInfo>();

        TableReader tr = new TableReader(NOVEL_KMER_MAP);

        for (Map<String, String> te : tr) {
            VariantInfo vi = allVariants.get(te.get("variantId"));

            CortexKmer kmer = new CortexKmer(te.get("kmer"));

            vis.put(kmer, vi);
        }

        return vis;
    }

    private String formatAttributes(Map<String, Object> attrs) {
        List<String> attrArray = new ArrayList<String>();

        for (String attr : attrs.keySet()) {
            String value = attrs.get(attr).toString();

            attrArray.add(attr + "=\"" + value + "\"");
        }

        return Joiner.on(" ").join(attrArray);
    }

    public void printGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, String prefix, boolean withText, boolean withPdf) {
        try {
            File f = new File(gout.getAbsolutePath() + "." + prefix + ".dot");
            File p = new File(gout.getAbsolutePath() + "." + prefix + ".pdf");

            PrintStream o = new PrintStream(f);

            String indent = "  ";

            o.println("digraph G {");

            if (withText) {
                o.println(indent + "node [ shape=box fontname=\"Courier New\" style=filled fillcolor=white ];");
            } else {
                o.println(indent + "node [ shape=point style=filled fillcolor=white ];");
            }

            for (AnnotatedVertex v : g.vertexSet()) {
                Map<String, Object> attrs = new TreeMap<String, Object>();

                if (!withText) {
                    attrs.put("label", "");
                }

                if (v.isNovel()) {
                    attrs.put("color", "red");
                    attrs.put("fillcolor", "red");
                    attrs.put("shape", "circle");
                }

                if (v.flagIsSet("start") || v.flagIsSet("end")) {
                    //attrs.put("label", v.getKmer());
                    attrs.put("fillcolor", v.flagIsSet("start") ? "orange" : "purple");
                    attrs.put("shape", "rect");
                }

                o.println(indent + "\"" + v.getKmer() + "\" [ " + formatAttributes(attrs) + " ];");
            }

            String[] colors = new String[] { "red", "blue", "green" };

            for (AnnotatedEdge e : g.edgeSet()) {
                String s = g.getEdgeSource(e).getKmer();
                String t = g.getEdgeTarget(e).getKmer();

                for (int c = 0; c < 3; c++) {
                    if (e.isPresent(c)) {
                        Map<String, Object> attrs = new TreeMap<String, Object>();
                        attrs.put("color", colors[c]);

                        o.println(indent + "\"" + s + "\" -> \"" + t + "\" [ " + formatAttributes(attrs) + " ];");
                    }
                }
            }

            o.println("}");

            o.close();

            if (withPdf) {
                //log.info("dot -Tpdf -o" + p.getAbsolutePath() + " " + f.getAbsolutePath());
                Runtime.getRuntime().exec("dot -Tpdf -o" + p.getAbsolutePath() + " " + f.getAbsolutePath());
            }
        } catch (FileNotFoundException e) {
            throw new IndianaException("File not found", e);
        } catch (IOException e) {
            throw new IndianaException("IO exception", e);
        }
    }

    private void addGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, DirectedGraph<String, DefaultEdge> g, int color, Map<CortexKmer, Boolean> novelKmers) {
        for (String v : g.vertexSet()) {
            AnnotatedVertex av = new AnnotatedVertex(v, novelKmers.containsKey(new CortexKmer(v)));

            a.addVertex(av);
        }

        for (DefaultEdge e : g.edgeSet()) {
            String s0 = g.getEdgeSource(e);
            String s1 = g.getEdgeTarget(e);

            CortexKmer ck0 = new CortexKmer(s0);
            CortexKmer ck1 = new CortexKmer(s1);

            AnnotatedVertex a0 = new AnnotatedVertex(s0, novelKmers.containsKey(ck0));
            AnnotatedVertex a1 = new AnnotatedVertex(s1, novelKmers.containsKey(ck1));

            //log.info("a0: {}", a0);
            //log.info("a1: {}", a1);

            if (!a.containsEdge(a0, a1)) {
                a.addEdge(a0, a1, new AnnotatedEdge());
            }

            a.getEdge(a0, a1).set(color, true);
        }
    }

    private Set<AnnotatedVertex> getListOfVerticesToRemove(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, DepthFirstIterator<AnnotatedVertex, AnnotatedEdge> dfs) {
        Set<AnnotatedVertex> toRemove = new HashSet<AnnotatedVertex>();

        while (dfs.hasNext()) {
            AnnotatedVertex v1 = dfs.next();

            if (!v1.isNovel() && !toRemove.contains(v1)) {
                Set<AnnotatedEdge> oes = a.outgoingEdgesOf(v1);

                for (AnnotatedEdge oe : oes) {
                    if (oe.isPresent(0) && (oe.isPresent(1) || oe.isPresent(2))) {
                        AnnotatedVertex vt = a.getEdgeTarget(oe);

                        DepthFirstIterator<AnnotatedVertex, AnnotatedEdge> subdfs = new DepthFirstIterator<AnnotatedVertex, AnnotatedEdge>(a, vt);

                        while (subdfs.hasNext()) {
                            AnnotatedVertex vs = subdfs.next();

                            toRemove.add(vs);
                        }
                    }
                }
            }
        }

        return toRemove;
    }

    private DirectedGraph<AnnotatedVertex, AnnotatedEdge> copyGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a) {
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> b = new DefaultDirectedGraph<AnnotatedVertex, AnnotatedEdge>(AnnotatedEdge.class);

        for (AnnotatedVertex av : a.vertexSet()) {
            AnnotatedVertex newAv = new AnnotatedVertex(av.getKmer(), av.isNovel());

            newAv.setFlags(av.getFlags());

            b.addVertex(newAv);
        }

        for (AnnotatedEdge ae : a.edgeSet()) {
            AnnotatedVertex vs = a.getEdgeSource(ae);
            AnnotatedVertex vt = a.getEdgeTarget(ae);

            b.addEdge(vs, vt, new AnnotatedEdge(ae.getPresence()));
        }

        return b;
    }

    private DirectedGraph<AnnotatedVertex, AnnotatedEdge> simplifyGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a) {
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> b = copyGraph(a);

        AnnotatedVertex thisVertex;
        Set<AnnotatedVertex> usedStarts = new HashSet<AnnotatedVertex>();

        Set<AnnotatedVertex> candidateStarts = new HashSet<AnnotatedVertex>();
        Set<AnnotatedVertex> candidateEnds = new HashSet<AnnotatedVertex>();
        for (AnnotatedVertex av : b.vertexSet()) {
            if (av.flagIsSet("start")) {
                candidateStarts.add(av);
            } else if (av.flagIsSet("end")) {
                candidateEnds.add(av);
            }
        }

        do {
            thisVertex = null;

            Iterator<AnnotatedVertex> vertices = b.vertexSet().iterator();

            while (thisVertex == null && vertices.hasNext()) {
                AnnotatedVertex tv = vertices.next();

                if (!usedStarts.contains(tv)) {
                    while (b.inDegreeOf(tv) == 1) {
                        AnnotatedVertex pv = b.getEdgeSource(b.incomingEdgesOf(tv).iterator().next());

                        if (b.outDegreeOf(pv) == 1 && pv.isNovel() == tv.isNovel()) {
                            tv = pv;
                        } else {
                            break;
                        }
                    }

                    thisVertex = tv;
                }
            }

            if (thisVertex != null) {
                usedStarts.add(thisVertex);

                AnnotatedVertex nextVertex = b.outDegreeOf(thisVertex) == 0 ? null : b.getEdgeTarget(b.outgoingEdgesOf(thisVertex).iterator().next());

                while (nextVertex != null && b.outDegreeOf(thisVertex) == 1 && b.inDegreeOf(nextVertex) == 1 && thisVertex.isNovel() == nextVertex.isNovel()) {
                    AnnotatedVertex sv = new AnnotatedVertex(thisVertex.getKmer() + nextVertex.getKmer().charAt(nextVertex.getKmer().length() - 1), thisVertex.isNovel());
                    //sv.setFlags(thisVertex.getFlags());

                    b.addVertex(sv);

                    Set<AnnotatedEdge> edgesToRemove = new HashSet<AnnotatedEdge>();

                    for (AnnotatedEdge e : b.incomingEdgesOf(thisVertex)) {
                        AnnotatedVertex pv = b.getEdgeSource(e);

                        b.addEdge(pv, sv, new AnnotatedEdge(e.getPresence()));

                        edgesToRemove.add(e);
                    }

                    for (AnnotatedEdge e : b.outgoingEdgesOf(nextVertex)) {
                        AnnotatedVertex nv = b.getEdgeTarget(e);

                        b.addEdge(sv, nv, new AnnotatedEdge(e.getPresence()));

                        edgesToRemove.add(e);
                    }

                    b.removeVertex(thisVertex);
                    b.removeVertex(nextVertex);
                    b.removeAllEdges(edgesToRemove);

                    thisVertex = sv;
                    nextVertex = b.outDegreeOf(thisVertex) == 0 ? null : b.getEdgeTarget(b.outgoingEdgesOf(thisVertex).iterator().next());
                }
            }
        } while (thisVertex != null);

        for (AnnotatedVertex av : b.vertexSet()) {
            for (AnnotatedVertex candidateStart : candidateStarts) {
                if (av.getKmer().contains(candidateStart.getKmer())) {
                    av.setFlag("start");
                }
            }

            for (AnnotatedVertex candidateEnd : candidateEnds) {
                if (av.getKmer().contains(candidateEnd.getKmer())) {
                    av.setFlag("end");
                }
            }
        }

        return b;
    }

    private DirectedGraph<AnnotatedVertex, AnnotatedEdge> removeOtherColors(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, int colorToRetain) {
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> b = copyGraph(a);

        Set<AnnotatedEdge> edgesToRemove = new HashSet<AnnotatedEdge>();

        for (AnnotatedEdge ae : b.edgeSet()) {
            for (int c = 0; c < 10; c++) {
                if (c != colorToRetain) {
                    ae.set(c, false);
                }
            }

            if (ae.isAbsent(colorToRetain)) {
                edgesToRemove.add(ae);
            }
        }

        b.removeAllEdges(edgesToRemove);

        Set<AnnotatedVertex> verticesToRemove = new HashSet<AnnotatedVertex>();
        for (AnnotatedVertex av : b.vertexSet()) {
            if (b.inDegreeOf(av) == 0 && b.outDegreeOf(av) == 0) {
                verticesToRemove.add(av);
            }
        }

        b.removeAllVertices(verticesToRemove);

        return b;
    }

    private String linearizePath(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, GraphPath<AnnotatedVertex, AnnotatedEdge> p) {
        StringBuilder sb = new StringBuilder(p.getStartVertex().getKmer());

        for (AnnotatedEdge ae : p.getEdgeList()) {
            AnnotatedVertex av = a.getEdgeTarget(ae);

            sb.append(av.getKmer().charAt(av.getKmer().length() - 1));
        }

        return sb.toString();
    }

    private DirectedGraph<AnnotatedVertex, AnnotatedEdge> trimGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, Set<AnnotatedVertex> candidateStarts, Set<AnnotatedVertex> candidateEnds) {
        /*
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> b = copyGraph(a);
        Set<AnnotatedVertex> verticesToRemove = new HashSet<AnnotatedVertex>();

        for (AnnotatedVertex candidateStart : candidateStarts) {
            Set<AnnotatedVertex> candidateVertices = new HashSet<AnnotatedVertex>(Graphs.predecessorListOf(b, candidateStart));

            boolean goodToRemove = true;
            for (AnnotatedVertex av : candidateVertices) {
                if (av.isNovel()) {
                    goodToRemove = false;
                }
            }

            if (goodToRemove) {
                verticesToRemove.addAll(candidateVertices);
            }
        }

        for (AnnotatedVertex candidateEnd : candidateEnds) {
            Set<AnnotatedVertex> candidateVertices = new HashSet<AnnotatedVertex>(Graphs.successorListOf(b, candidateEnd));

            boolean goodToRemove = true;
            for (AnnotatedVertex av : candidateVertices) {
                if (av.isNovel()) {
                    goodToRemove = false;
                }
            }

            if (goodToRemove) {
                verticesToRemove.addAll(candidateVertices);
            }
        }

        Set<AnnotatedEdge> edgesToRemove = new HashSet<AnnotatedEdge>();

        for (AnnotatedVertex av : verticesToRemove) {
            for (AnnotatedEdge ae : b.outgoingEdgesOf(av)) {
                AnnotatedVertex nav = b.getEdgeTarget(ae);

                if (verticesToRemove.contains(nav)) {
                    edgesToRemove.add(ae);
                }
            }

            for (AnnotatedEdge ae : b.incomingEdgesOf(av)) {
                AnnotatedVertex pav = b.getEdgeSource(ae);

                if (verticesToRemove.contains(pav)) {
                    edgesToRemove.add(ae);
                }
            }
        }

        b.removeAllVertices(verticesToRemove);
        b.removeAllEdges(edgesToRemove);

        return b;
        */

        return a;
    }

    private Pair<String, String> computeMinWeightPaths(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, int color, String stretch, Map<CortexKmer, Boolean> novelKmers) {
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> b = copyGraph(a);

        log.info("    a={} {}, b={} {}", a.vertexSet().size(), a.edgeSet().size(), b.vertexSet().size(), b.edgeSet().size());

        AnnotateStartsAndEnds annotateStartsAndEnds = new AnnotateStartsAndEnds(color, stretch, novelKmers, b).invoke();
        Set<AnnotatedVertex> candidateStarts = annotateStartsAndEnds.getCandidateStarts();
        Set<AnnotatedVertex> candidateEnds = annotateStartsAndEnds.getCandidateEnds();

        DirectedGraph<AnnotatedVertex, AnnotatedEdge> b0 = removeOtherColors(b, 0);
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> bc = removeOtherColors(b, color);

        //printGraph(simplifyGraph(trimGraph(b, candidateStarts, candidateEnds)), "monocolor.all",      false, true);
        //printGraph(simplifyGraph(trimGraph(b0, candidateStarts, candidateEnds)), "monocolor." + 0,     false, true);
        //printGraph(simplifyGraph(trimGraph(bc, candidateStarts, candidateEnds)), "monocolor." + color, false, true);

        GraphPath<AnnotatedVertex, AnnotatedEdge> p0 = computeMinWeightPath(b0, candidateEnds, candidateStarts);
        GraphPath<AnnotatedVertex, AnnotatedEdge> pc = computeMinWeightPath(bc, candidateEnds, candidateStarts);

        return new Pair<String, String>(
                p0 == null ? "" : linearizePath(b0, p0),
                pc == null ? "" : linearizePath(bc, pc)
        );
    }

    private GraphPath<AnnotatedVertex, AnnotatedEdge> computeMinWeightPath(DirectedGraph<AnnotatedVertex, AnnotatedEdge> b, Set<AnnotatedVertex> candidateEnds, Set<AnnotatedVertex> candidateStarts) {
        double minPathLength = Double.MAX_VALUE;
        GraphPath<AnnotatedVertex, AnnotatedEdge> minWeightPath = null;

        for (AnnotatedVertex sv : candidateStarts) {
            for (AnnotatedVertex ev : candidateEnds) {
                DijkstraShortestPath<AnnotatedVertex, AnnotatedEdge> dsp = new DijkstraShortestPath<AnnotatedVertex, AnnotatedEdge>(b, sv, ev);

                GraphPath<AnnotatedVertex, AnnotatedEdge> path = dsp.getPath();
                double pathLength = dsp.getPathLength();

                if (pathLength < minPathLength) {
                    minPathLength = pathLength;
                    minWeightPath = path;
                }
            }
        }

        return minWeightPath;
    }

    private boolean isConnected(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, int color, String stretch, Map<CortexKmer, Boolean> novelKmers) {
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> b = copyGraph(a);

        AnnotateStartsAndEnds annotateStartsAndEnds = new AnnotateStartsAndEnds(color, stretch, novelKmers, b).invoke();
        Set<AnnotatedVertex> candidateStarts = annotateStartsAndEnds.getCandidateStarts();
        Set<AnnotatedVertex> candidateEnds = annotateStartsAndEnds.getCandidateEnds();

        b = removeOtherColors(b, color);

        ConnectivityInspector<AnnotatedVertex, AnnotatedEdge> ci = new ConnectivityInspector<AnnotatedVertex, AnnotatedEdge>(b);
        //StrongConnectivityInspector<AnnotatedVertex, AnnotatedEdge> ci = new StrongConnectivityInspector<AnnotatedVertex, AnnotatedEdge>(b);

        for (AnnotatedVertex candidateStart : candidateStarts) {
            for (AnnotatedVertex candidateEnd : candidateEnds) {
                if (ci.pathExists(candidateStart, candidateEnd)) {
                    log.info("s: {}, e: {}", candidateStart, candidateEnd);

                    return true;
                }
            }
        }

        return false;
    }

    private Pair<String, String> getAlleles(Pair<String, String> minWeightPaths) {
        int s, e0 = minWeightPaths.getFirst().length() - 1, e1 = minWeightPaths.getSecond().length() - 1;

        for (s = 0; s < (minWeightPaths.getFirst().length() < minWeightPaths.getSecond().length() ? minWeightPaths.getFirst().length() : minWeightPaths.getSecond().length()) && minWeightPaths.getFirst().charAt(s) == minWeightPaths.getSecond().charAt(s); s++) {}

        while (e0 > s && e1 > s && minWeightPaths.getFirst().charAt(e0) == minWeightPaths.getSecond().charAt(e1)) {
            e0--;
            e1--;
        }

        return new Pair<String, String>(minWeightPaths.getFirst().substring(s, e0 + 1), minWeightPaths.getSecond().substring(s, e1 + 1));
    }

    private VariantContext getVariant(Pair<String, String> minWeightPaths, String stretch) {
        int s, e0 = minWeightPaths.getFirst().length() - 1, e1 = minWeightPaths.getSecond().length() - 1;

        for (s = 0; s < (minWeightPaths.getFirst().length() < minWeightPaths.getSecond().length() ? minWeightPaths.getFirst().length() : minWeightPaths.getSecond().length()) && minWeightPaths.getFirst().charAt(s) == minWeightPaths.getSecond().charAt(s); s++) {}

        while (e0 > s && e1 > s && minWeightPaths.getFirst().charAt(e0) == minWeightPaths.getSecond().charAt(e1)) {
            e0--;
            e1--;
        }

        VariantContext vc = (new VariantContextBuilder())
                .chr("unknown")
                .start(s)
                .stop(e1)
                .alleles(minWeightPaths.getSecond().substring(s, e1 + 1), minWeightPaths.getFirst().substring(s, e0 + 1))
                .attribute("NOVEL_STRETCH", stretch)
                .attribute("CHILD_STRETCH", minWeightPaths.getFirst())
                .attribute("PARENT_STRETCH", minWeightPaths.getSecond())
                .make();

        return vc;
    }

    @Override
    public void execute() {
        TableWriter tw = new TableWriter(out);

        Map<CortexKmer, VariantInfo> vis = loadNovelKmerMap();

        Map<CortexKmer, Boolean> novelKmers = new HashMap<CortexKmer, Boolean>();
        for (CortexRecord cr : NOVEL) {
            novelKmers.put(new CortexKmer(cr.getKmerAsString()), true);
        }

        int totalNovelKmersUsed = 0;
        int stretchNum = 0;
        for (CortexKmer novelKmer : novelKmers.keySet()) {
            if (novelKmers.get(novelKmer)) {
                int novelKmersUsed = 0;

                // first, explore each color and bring the local subgraphs into memory
                DirectedGraph<String, DefaultEdge> sg0 = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);
                DirectedGraph<String, DefaultEdge> sg1 = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);
                DirectedGraph<String, DefaultEdge> sg2 = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);

                String stretch = CortexUtils.getSeededStretch(GRAPH, novelKmer.getKmerAsString(), 0);

                log.info("  stretch: {} ({})", stretchNum, stretch.length());
                log.info("    processing child graph...");

                for (int i = 0; i <= stretch.length() - novelKmer.length(); i++) {
                    String kmer = stretch.substring(i, i + novelKmer.length());

                    if (!sg0.containsVertex(kmer)) {
                        Graphs.addGraph(sg0, CortexUtils.dfs(GRAPH, kmer, 0, null, new AbstractTraversalStopper() {
                            @Override
                            public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<String, DefaultEdge> g, int depth) {
                                return false;
                            }

                            @Override
                            public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<String, DefaultEdge> g, int junctions) {
                                return junctions >= maxJunctions;
                            }
                        }));
                    }
                }

                for (int c = 1; c <= 2; c++) {
                    log.info("    processing parent {} graph... ({} kmers)", c, sg0.vertexSet().size());

                    for (String kmer : sg0.vertexSet()) {
                        DirectedGraph<String, DefaultEdge> sg = (c == 1) ? sg1 : sg2;

                        if (!sg.containsVertex(kmer)) {
                            TraversalStopper stopper = new AbstractTraversalStopper() {
                                @Override
                                public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<String, DefaultEdge> g, int junctions) {
                                    String fw = cr.getKmerAsString();
                                    String rc = SequenceUtils.reverseComplement(fw);

                                    if (g.containsVertex(fw) || g.containsVertex(rc)) {
                                        if (junctions < distanceToGoal) {
                                            distanceToGoal = junctions;
                                        }

                                        return true;
                                    }

                                    return false;
                                }

                                @Override
                                public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<String, DefaultEdge> g, int junctions) {
                                    return junctions >= maxJunctions;
                                }
                            };

                            Graphs.addGraph(sg, CortexUtils.dfs(GRAPH, kmer, c, sg0, stopper));
                        }
                    }
                }

                log.info("    combining graphs...");

                // Now, combine them all into an annotated graph
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> ag = new DefaultDirectedGraph<AnnotatedVertex, AnnotatedEdge>(AnnotatedEdge.class);

                addGraph(ag, sg0, 0, novelKmers);
                addGraph(ag, sg1, 1, novelKmers);
                addGraph(ag, sg2, 2, novelKmers);

                printGraph(simplifyGraph(ag), "debug", false, false);

                Set<VariantInfo> relevantVis = new HashSet<VariantInfo>();
                for (AnnotatedVertex ak : ag.vertexSet()) {
                    CortexKmer ck = new CortexKmer(ak.getKmer());

                    if (ak.isNovel() && vis.containsKey(ck)) {
                        VariantInfo vi = vis.get(novelKmer);

                        relevantVis.add(vi);
                    }
                }

                // Perform novel kmer accounting
                for (AnnotatedVertex ak : ag.vertexSet()) {
                    CortexKmer ck = new CortexKmer(ak.getKmer());

                    if (novelKmers.containsKey(ck) && novelKmers.get(ck)) {
                        totalNovelKmersUsed++;
                        novelKmersUsed++;
                        novelKmers.put(ck, false);
                    }
                }

                log.info("    novelKmers: {}, totalNovelKmersUsed: {}, totalNovelKmers: {}", novelKmersUsed, totalNovelKmersUsed, novelKmers.size());

                // Now, find min-weight paths for each color vs child
                boolean c1 = isConnected(ag, 1, stretch, novelKmers);
                Pair<String, String> p1 = computeMinWeightPaths(ag, 1, stretch, novelKmers);
                VariantContext vc1 = (new VariantContextBuilder())
                        .chr("unknown")
                        .start(1)
                        .stop(1)
                        .alleles("A", "N")
                        .filter("TRAVERSAL_FAILED")
                        .attribute("NOVEL_STRETCH", stretch)
                        .make();
                VariantContext vc2 = (new VariantContextBuilder())
                        .chr("unknown")
                        .start(1)
                        .stop(1)
                        .alleles("A", "N")
                        .filter("TRAVERSAL_FAILED")
                        .attribute("NOVEL_STRETCH", stretch)
                        .make();

                boolean c2 = isConnected(ag, 2, stretch, novelKmers);
                Pair<String, String> p2 = computeMinWeightPaths(ag, 2, stretch, novelKmers);

                log.info("    vis: {}", relevantVis);
                log.info("     c1: {}", c1);
                log.info("     n1: {} ({})", p1.getFirst(), p1.getFirst().length());
                log.info("     p1: {} ({})", p1.getSecond(), p1.getSecond().length());

                if (p1.getFirst() != null && p1.getSecond() != null && p1.getFirst().length() > 0 && p1.getSecond().length() > 0 && p1.getFirst().charAt(0) == p1.getSecond().charAt(0) && !p1.getFirst().equals(p1.getSecond())) {
                    Pair<String, String> a1 = getAlleles(p1);
                    vc1 = getVariant(p1, stretch);

                    log.info("     a1: {}", a1.getFirst());
                    log.info("     a1: {}", a1.getSecond());
                    log.info("     v1: {}", vc1);
                }

                log.info("     c2: {}", c2);
                log.info("     n2: {} ({})", p2.getFirst(), p2.getFirst().length());
                log.info("     p2: {} ({})", p2.getSecond(), p2.getSecond().length());

                if (p2.getFirst() != null && p2.getSecond() != null && p2.getFirst().length() > 0 && p2.getSecond().length() > 0 && p2.getFirst().charAt(0) == p2.getSecond().charAt(0) && !p2.getFirst().equals(p2.getSecond())) {
                    Pair<String, String> a2 = getAlleles(p2);
                    vc2 = getVariant(p2, stretch);

                    log.info("     a2: {}", a2.getFirst());
                    log.info("     a2: {}", a2.getSecond());
                    log.info("     v2: {}", vc2);
                }

                Map<String, String> te = new LinkedHashMap<String, String>();

                te.put("variantId", "unknown");
                te.put("knownType", "unknown");
                te.put("knownRef", "unknown");
                te.put("knownAlt", "unknown");

                for (VariantInfo vi : relevantVis) {
                    te.put("variantId", vi.variantId);
                    te.put("knownType", vi.type);
                    te.put("knownRef", vi.ref);
                    te.put("knownAlt", vi.alt == null ? "." : vi.alt);
                }

                if (vc1 != null) {
                    String ref = vc1.getReference().getBaseString();
                    String alt = vc1.getAlternateAllele(0).getBaseString();

                    if (!te.get("knownRef").equals("unknown") && te.get("knownRef").length() == vc1.getReference().length() && !te.get("knownRef").equals(vc1.getReference().getBaseString()) && !te.get("knownAlt").equals(vc1.getAlternateAllele(0).getBaseString())) {
                        int pos = vc1.getStart();
                        int refLength = vc1.getReference().length();
                        int altLength = vc1.getAlternateAllele(0).length();
                        boolean found = false;

                        while (pos >= 0 && pos + refLength < p1.getFirst().length() && pos + altLength < p1.getSecond().length()) {
                            String refFw = p1.getSecond().substring(pos, pos + refLength);
                            String refRc = SequenceUtils.reverseComplement(refFw);

                            String altFw = p1.getFirst().substring(pos, pos + altLength);
                            String altRc = SequenceUtils.reverseComplement(altFw);

                            if (te.get("knownRef").equals(refFw) && te.get("knownAlt").equals(altFw)) {
                                ref = refFw;
                                alt = altFw;
                                found = true;
                                break;
                            } else if (te.get("knownRef").equals(refRc) && te.get("knownAlt").equals(altRc)) {
                                ref = refRc;
                                alt = altRc;
                                found = true;
                                break;
                            }

                            pos--;
                        }

                        if (!found) {
                            pos = vc1.getStart();

                            while (pos >= 0 && pos + refLength < p1.getFirst().length() && pos + altLength < p1.getSecond().length()) {
                                String refFw = p1.getSecond().substring(pos, pos + refLength);
                                String refRc = SequenceUtils.reverseComplement(refFw);

                                String altFw = p1.getFirst().substring(pos, pos + altLength);
                                String altRc = SequenceUtils.reverseComplement(altFw);

                                if (te.get("knownRef").equals(refFw) && te.get("knownAlt").equals(altFw)) {
                                    ref = refFw;
                                    alt = altFw;
                                    break;
                                } else if (te.get("knownRef").equals(refRc) && te.get("knownAlt").equals(altRc)) {
                                    ref = refRc;
                                    alt = altRc;
                                    break;
                                }

                                pos++;
                            }
                        }
                    }

                    te.put("v1Type", vc1.getType().name());
                    te.put("v1Ref", ref);
                    te.put("v1Alt", alt);
                    te.put("v1Filter", String.valueOf(vc1.isFiltered()));
                } else {
                    te.put("v1Type", "NA");
                    te.put("v1Ref", "NA");
                    te.put("v1Alt", "NA");
                    te.put("v1Filter", "NA");
                }

                if (vc2 != null) {
                    String ref = vc2.getReference().getBaseString();
                    String alt = vc2.getAlternateAllele(0).getBaseString();

                    if (!te.get("knownRef").equals("unknown") && te.get("knownRef").length() == vc2.getReference().length() && !te.get("knownRef").equals(vc2.getReference().getBaseString()) && !te.get("knownAlt").equals(vc2.getAlternateAllele(0).getBaseString())) {
                        int pos = vc2.getStart();
                        int refLength = vc2.getReference().length();
                        int altLength = vc2.getAlternateAllele(0).length();
                        boolean found = false;

                        while (pos >= 0 && pos + refLength < p2.getFirst().length() && pos + altLength < p2.getSecond().length()) {
                            String refFw = p2.getSecond().substring(pos, pos + vc2.getReference().length());
                            String refRc = SequenceUtils.reverseComplement(refFw);

                            String altFw = p2.getFirst().substring(pos, pos + vc2.getAlternateAllele(0).length());
                            String altRc = SequenceUtils.reverseComplement(altFw);

                            if (te.get("knownRef").equals(refFw) && te.get("knownAlt").equals(altFw)) {
                                ref = refFw;
                                alt = altFw;
                                found = true;
                                break;
                            } else if (te.get("knownRef").equals(refRc) && te.get("knownAlt").equals(altRc)) {
                                ref = refRc;
                                alt = altRc;
                                found = true;
                                break;
                            }

                            pos--;
                        }

                        if (!found) {
                            pos = vc2.getStart();

                            while (pos >= 0 && pos + refLength < p2.getFirst().length() && pos + altLength < p2.getSecond().length()) {
                                String refFw = p2.getSecond().substring(pos, pos + refLength);
                                String refRc = SequenceUtils.reverseComplement(refFw);

                                String altFw = p2.getFirst().substring(pos, pos + altLength);
                                String altRc = SequenceUtils.reverseComplement(altFw);

                                if (te.get("knownRef").equals(refFw) && te.get("knownAlt").equals(altFw)) {
                                    ref = refFw;
                                    alt = altFw;
                                    break;
                                } else if (te.get("knownRef").equals(refRc) && te.get("knownAlt").equals(altRc)) {
                                    ref = refRc;
                                    alt = altRc;
                                    break;
                                }

                                pos++;
                            }
                        }
                    }

                    te.put("v2Type", vc2.getType().name());
                    te.put("v2Ref", ref);
                    te.put("v2Alt", alt);
                    te.put("v2Filter", String.valueOf(vc2.isFiltered()));
                } else {
                    te.put("v2Type", "NA");
                    te.put("v2Ref", "NA");
                    te.put("v2Alt", "NA");
                    te.put("v2Filter", "NA");
                }

                tw.addEntry(te);
                out.flush();

                stretchNum++;
            }
        }

        log.info("Num stretches: {}", stretchNum);
    }

    private class AnnotateStartsAndEnds {
        private int color;
        private String stretch;
        private Map<CortexKmer, Boolean> novelKmers;
        private DirectedGraph<AnnotatedVertex, AnnotatedEdge> b;
        private Set<AnnotatedVertex> candidateEnds;
        private Set<AnnotatedVertex> candidateStarts;

        public AnnotateStartsAndEnds(int color, String stretch, Map<CortexKmer, Boolean> novelKmers, DirectedGraph<AnnotatedVertex, AnnotatedEdge> b) {
            this.color = color;
            this.stretch = stretch;
            this.novelKmers = novelKmers;
            this.b = b;
        }

        public Set<AnnotatedVertex> getCandidateEnds() {
            return candidateEnds;
        }

        public Set<AnnotatedVertex> getCandidateStarts() {
            return candidateStarts;
        }

        public AnnotateStartsAndEnds invoke() {
            int kmerLength = novelKmers.keySet().iterator().next().length();

            String fk = stretch.substring(1, kmerLength + 1);
            AnnotatedVertex afk = new AnnotatedVertex(fk, novelKmers.containsKey(new CortexKmer(fk)));

            String lk = stretch.substring(stretch.length() - kmerLength - 1, stretch.length() - 1);
            AnnotatedVertex alk = new AnnotatedVertex(lk, novelKmers.containsKey(new CortexKmer(lk)));

            candidateEnds = new HashSet<AnnotatedVertex>();
            DepthFirstIterator<AnnotatedVertex, AnnotatedEdge> dfsLk = new DepthFirstIterator<AnnotatedVertex, AnnotatedEdge>(b, afk);
            while (dfsLk.hasNext()) {
                AnnotatedVertex av = dfsLk.next();

                if (b.inDegreeOf(av) > 1) {
                    Set<AnnotatedEdge> aes = b.incomingEdgesOf(av);

                    boolean adjacentToNovelty = false;
                    AnnotatedEdge edgeToNovelty = null;

                    for (AnnotatedEdge ae : aes) {
                        if (ae.isPresent(0) && ae.isAbsent(1) && ae.isAbsent(2)) {
                            adjacentToNovelty = true;
                            edgeToNovelty = ae;
                        }
                    }

                    if (adjacentToNovelty) {
                        for (AnnotatedEdge ae : aes) {
                            if (!ae.equals(edgeToNovelty)) {
                                if (ae.isPresent(color)) {
                                    candidateEnds.add(av);

                                    //av.setFlag("end");

                                    //
                                    for (AnnotatedVertex ava : b.vertexSet()) {
                                        if (av.equals(ava)) {
                                            ava.setFlag("end");
                                        }
                                    }
                                    //
                                }
                            }
                        }
                    }
                }
            }

            candidateStarts = new HashSet<AnnotatedVertex>();
            DirectedGraph<AnnotatedVertex, AnnotatedEdge> rag = new EdgeReversedGraph<AnnotatedVertex, AnnotatedEdge>(b);
            DepthFirstIterator<AnnotatedVertex, AnnotatedEdge> dfsFk = new DepthFirstIterator<AnnotatedVertex, AnnotatedEdge>(rag, alk);
            while (dfsFk.hasNext()) {
                AnnotatedVertex av = dfsFk.next();

                if (rag.inDegreeOf(av) > 1) {
                    Set<AnnotatedEdge> aes = rag.incomingEdgesOf(av);

                    boolean adjacentToNovelty = false;
                    AnnotatedEdge edgeToNovelty = null;

                    for (AnnotatedEdge ae : aes) {
                        if (ae.isPresent(0) && ae.isAbsent(1) && ae.isAbsent(2)) {
                            adjacentToNovelty = true;
                            edgeToNovelty = ae;
                        }
                    }

                    if (adjacentToNovelty) {
                        for (AnnotatedEdge ae : aes) {
                            if (!ae.equals(edgeToNovelty)) {
                                if (ae.isPresent(color)) {
                                    candidateStarts.add(av);

                                    //av.setFlag("start");

                                    //
                                    for (AnnotatedVertex ava : b.vertexSet()) {
                                        if (av.equals(ava)) {
                                            ava.setFlag("start");
                                        }
                                    }
                                    //
                                }
                            }
                        }
                    }
                }
            }
            return this;
        }
    }
}
