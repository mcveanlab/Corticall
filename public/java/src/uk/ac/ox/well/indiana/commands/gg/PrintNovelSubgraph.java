package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
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
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

public class PrintNovelSubgraph extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="graphRaw", shortName="gr", doc="Graph (raw)", required=false)
    public CortexGraph GRAPH_RAW;

    @Argument(fullName="novelGraph", shortName="n", doc="Graph of novel kmers")
    public CortexGraph NOVEL;

    @Argument(fullName="novelKmerMap", shortName="m", doc="Novel kmer map")
    public File NOVEL_KMER_MAP;

    @Output
    public File out;

    private class VariantInfo {
        public String variantId;
        public String vclass;
        public String vchr;
        public int vstart;
        public int vstop;

        @Override
        public String toString() {
            return "VariantInfo{" +
                    "variantId='" + variantId + '\'' +
                    ", vclass='" + vclass + '\'' +
                    ", vchr='" + vchr + '\'' +
                    ", vstart=" + vstart +
                    ", vstop=" + vstop +
                    '}';
        }
    }

    private Map<CortexKmer, VariantInfo> loadNovelKmerMap() {
        TableReader tr = new TableReader(NOVEL_KMER_MAP);

        Map<CortexKmer, VariantInfo> vis = new HashMap<CortexKmer, VariantInfo>();

        for (Map<String, String> te : tr) {
            VariantInfo vi = new VariantInfo();

            if (!te.get("vstart").equals("NA")) {
                vi.variantId = te.get("variantId");
                vi.vclass = te.get("vclass");
                vi.vchr = te.get("vchr");
                vi.vstart = Integer.valueOf(te.get("vstart"));
                vi.vstop = Integer.valueOf(te.get("vstop"));

                CortexKmer kmer = new CortexKmer(te.get("kmer"));

                vis.put(kmer, vi);
            }
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

    public void printGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, String prefix) {
        try {
            File f = new File(out.getAbsolutePath() + "." + prefix + ".dot");
            File p = new File(out.getAbsolutePath() + "." + prefix + ".png");

            PrintStream o = new PrintStream(f);

            String indent = "  ";

            o.println("digraph G {");

            for (AnnotatedVertex v : g.vertexSet()) {
                Map<String, Object> attrs = new TreeMap<String, Object>();
                //attrs.put("label", "");
                attrs.put("fillcolor", v.isNovel() ? "red" : "white");
                attrs.put("style", "filled");

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

            log.info("dot -Tpng -o" + p.getAbsolutePath() + " " + f.getAbsolutePath());
            Runtime.getRuntime().exec("dot -Tpng -o" + p.getAbsolutePath() + " " + f.getAbsolutePath());
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
    private void simplifyGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a) {
        AnnotatedVertex thisVertex;
        Set<AnnotatedVertex> usedStarts = new HashSet<AnnotatedVertex>();

        do {
            thisVertex = null;

            Iterator<AnnotatedVertex> vertices = a.vertexSet().iterator();

            while (thisVertex == null && vertices.hasNext()) {
                AnnotatedVertex tv = vertices.next();

                if (!usedStarts.contains(tv)) {
                    while (a.inDegreeOf(tv) == 1) {
                        AnnotatedVertex pv = a.getEdgeSource(a.incomingEdgesOf(tv).iterator().next());

                        if (a.outDegreeOf(pv) == 1 && pv.isNovel() == tv.isNovel()) {
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

                AnnotatedVertex nextVertex = a.outDegreeOf(thisVertex) == 0 ? null : a.getEdgeTarget(a.outgoingEdgesOf(thisVertex).iterator().next());

                /*
                if (nextVertex == null || a.outDegreeOf(thisVertex) != 1 || a.inDegreeOf(nextVertex) != 1 || thisVertex.isNovel() != nextVertex.isNovel()) {
                    //thisVertex = null;
                    break;
                }
                */

                while (nextVertex != null && a.outDegreeOf(thisVertex) == 1 && a.inDegreeOf(nextVertex) == 1 && thisVertex.isNovel() == nextVertex.isNovel()) {
                    AnnotatedVertex sv = new AnnotatedVertex(thisVertex.getKmer() + nextVertex.getKmer().charAt(nextVertex.getKmer().length() - 1), thisVertex.isNovel());

                    a.addVertex(sv);

                    Set<AnnotatedEdge> edgesToRemove = new HashSet<AnnotatedEdge>();

                    for (AnnotatedEdge e : a.incomingEdgesOf(thisVertex)) {
                        AnnotatedVertex pv = a.getEdgeSource(e);

                        a.addEdge(pv, sv, new AnnotatedEdge(e.getPresence()));

                        edgesToRemove.add(e);
                    }

                    for (AnnotatedEdge e : a.outgoingEdgesOf(nextVertex)) {
                        AnnotatedVertex nv = a.getEdgeTarget(e);

                        a.addEdge(sv, nv, new AnnotatedEdge(e.getPresence()));

                        edgesToRemove.add(e);
                    }

                    a.removeVertex(thisVertex);
                    a.removeVertex(nextVertex);
                    a.removeAllEdges(edgesToRemove);

                    thisVertex = sv;
                    nextVertex = a.outDegreeOf(thisVertex) == 0 ? null : a.getEdgeTarget(a.outgoingEdgesOf(thisVertex).iterator().next());

                    log.info("sv: {}", sv);
                }

                log.info("Loop done");
            }
        } while (thisVertex != null);
    }

    /*
    private void simplifyGraph2(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a) {
        AnnotatedVertex startVertex;

        do {
            startVertex = null;
            for (AnnotatedVertex av : a.vertexSet()) {
                //if (a.inDegreeOf(av) != 1 && a.outDegreeOf(av) == 1 && a.getEdgeTarget(a.outgoingEdgesOf(av).iterator().next()).isNovel() == av.isNovel()) {
                boolean prevKmersAreNovel = false, nextKmersAreNovel = false;

                for (AnnotatedEdge e : a.incomingEdgesOf(av)) {
                    AnnotatedVertex pv = a.getEdgeSource(e);

                    if (pv.isNovel()) { prevKmersAreNovel = true; }
                }

                for (AnnotatedEdge e : a.outgoingEdgesOf(av)) {
                    AnnotatedVertex nv = a.getEdgeTarget(e);

                    if (nv.isNovel()) { nextKmersAreNovel = true; }
                }

                boolean startOfNovelty = !prevKmersAreNovel && nextKmersAreNovel;

                //if (a.outDegreeOf(av) == 1 && (a.inDegreeOf(av) != 1 || startOfNovelty)) {
                if (a.outDegreeOf(av) == 1 && a.inDegreeOf(av) != 1) {
                    startVertex = av;
                    break;
                }
            }

            if (startVertex != null && a.containsVertex(startVertex)) {
                Set<AnnotatedVertex> verticesBeingProcessed = new HashSet<AnnotatedVertex>();
                Set<AnnotatedEdge> edgesBeingProcessed = new HashSet<AnnotatedEdge>();

                StringBuilder supernode = new StringBuilder(startVertex.getKmer());

                AnnotatedVertex v0 = startVertex;
                Set<AnnotatedEdge> eps = a.incomingEdgesOf(v0);

                while (a.outDegreeOf(v0) == 1 && v0.isNovel() == startVertex.isNovel()) {
                    AnnotatedEdge ae = a.outgoingEdgesOf(v0).iterator().next();
                    AnnotatedVertex vn = a.getEdgeTarget(ae);

                    supernode.append(vn.getKmer().charAt(vn.getKmer().length() - 1));

                    verticesBeingProcessed.add(v0);
                    edgesBeingProcessed.add(ae);

                    v0 = vn;
                }

                Set<AnnotatedEdge> ens = a.outgoingEdgesOf(v0);

                AnnotatedVertex sv = new AnnotatedVertex(supernode.toString());
                if (startVertex.isNovel()) {
                    sv.setNovel();
                }

                a.addVertex(sv);

                for (AnnotatedEdge e : eps) {
                    AnnotatedVertex vs = a.getEdgeSource(e);
                    //a.addEdge(vs, sv, e);

                    a.addEdge(vs, sv, new AnnotatedEdge(e.getPresence()));
                }

                for (AnnotatedEdge e : ens) {
                    AnnotatedVertex vt = a.getEdgeTarget(e);
                    //a.addEdge(sv, vt, e);

                    a.addEdge(sv, vt, new AnnotatedEdge(e.getPresence()));
                }

                a.removeAllVertices(verticesBeingProcessed);
                a.removeAllEdges(edgesBeingProcessed);

                //log.info("  {} {}", supernode.length(), verticesBeingProcessed.size());
            }
        } while (startVertex != null);
    }
    */

    private void trimGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> af) {
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> ar = new EdgeReversedGraph<AnnotatedVertex, AnnotatedEdge>(af);

        Set<AnnotatedVertex> toRemove = new HashSet<AnnotatedVertex>();

        for (AnnotatedVertex av : af.vertexSet()) {
            if (av.isNovel()) {
                AnnotatedVertex v0 = av;

                DepthFirstIterator<AnnotatedVertex, AnnotatedEdge> dfsf = new DepthFirstIterator<AnnotatedVertex, AnnotatedEdge>(af, v0);
                DepthFirstIterator<AnnotatedVertex, AnnotatedEdge> dfsr = new DepthFirstIterator<AnnotatedVertex, AnnotatedEdge>(ar, v0);

                Set<AnnotatedVertex> toRemoveFwd = getListOfVerticesToRemove(af, dfsf);
                Set<AnnotatedVertex> toRemoveRev = getListOfVerticesToRemove(ar, dfsr);

                toRemoveFwd.remove(v0);
                toRemoveRev.remove(v0);

                //log.info("  fwd: {}/{}", toRemoveFwd.size(), af.vertexSet().size());
                //log.info("  rev: {}/{}", toRemoveRev.size(), af.vertexSet().size());

                toRemove.addAll(toRemoveFwd);
                toRemove.addAll(toRemoveRev);
            }
        }

        af.removeAllVertices(toRemove);
    }

    @Override
    public void execute() {
        Map<CortexKmer, VariantInfo> vis = loadNovelKmerMap();

        Map<CortexKmer, Boolean> novelKmers = new HashMap<CortexKmer, Boolean>();
        for (CortexRecord cr : NOVEL) {
            novelKmers.put(new CortexKmer(cr.getKmerAsString()), true);
        }

        int totalNovelKmersUsed = 0;
        int stretchNum = 0;
        for (CortexKmer novelKmer : novelKmers.keySet()) {
            if (novelKmers.get(novelKmer)) {
                VariantInfo vi = vis.get(novelKmer);

                if (vi != null && vi.vclass.equals("SNP")) {
                    log.info("vi: {}", vi);

                    int novelKmersUsed = 0;

                    // do stuff
                    DirectedGraph<String, DefaultEdge> sg0 = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);
                    DirectedGraph<String, DefaultEdge> sg1 = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);
                    DirectedGraph<String, DefaultEdge> sg2 = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);

                    String stretch = CortexUtils.getSeededStretch(GRAPH, novelKmer.getKmerAsString(), 0);

                    log.info("  stretch: {} ({})", stretchNum, stretch.length());
                    log.info("    processing child graph...");

                    /*
                    for (int i = 0; i <= stretch.length() - novelKmer.length(); i++) {
                        String kmer = stretch.substring(i, i + novelKmer.length());

                        if (!sg0.containsVertex(kmer)) {
                            Graphs.addGraph(sg0, CortexUtils.dfs(GRAPH, kmer, 0, null, new AbstractTraversalStopper() {
                                @Override
                                public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<String, DefaultEdge> g, int depth) {
                                    //return cr.getCoverage(1) != 0 || cr.getCoverage(2) != 0;
                                    return false;
                                }

                                @Override
                                public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<String, DefaultEdge> g, int junctions) {
                                    return junctions >= maxJunctions;
                                }
                            }));
                        }
                    }
                    */

                    for (int i = 0; i <= stretch.length() - novelKmer.length() - 1; i++) {
                        String k0 = stretch.substring(i, i + novelKmer.length());
                        String k1 = stretch.substring(i + 1, i + 1 + novelKmer.length());

                        sg0.addVertex(k0);
                        sg0.addVertex(k1);
                        sg0.addEdge(k0, k1);
                    }

                    //
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
                    //

                    log.info("    combining graphs...");

                    DirectedGraph<AnnotatedVertex, AnnotatedEdge> ag = new DefaultDirectedGraph<AnnotatedVertex, AnnotatedEdge>(AnnotatedEdge.class);

                    addGraph(ag, sg0, 0, novelKmers);
                    addGraph(ag, sg1, 1, novelKmers);
                    addGraph(ag, sg2, 2, novelKmers);

                    String prefix = String.format("%s.%s.%s.%d-%d", vi.variantId, vi.vclass, vi.vchr, vi.vstart, vi.vstop);

                    //printGraph(ag, prefix);


                    //log.info("    trimming graph...");
                    //trimGraph(ag);

                    printGraph(ag, prefix);

                    log.info("    simplifying graph...");
                    simplifyGraph(ag);

                    printGraph(ag, prefix + ".simplified");

                    for (AnnotatedVertex ak : ag.vertexSet()) {
                        CortexKmer ck = new CortexKmer(ak.getKmer());

                        if (novelKmers.containsKey(ck) && novelKmers.get(ck)) {
                            totalNovelKmersUsed++;
                            novelKmersUsed++;
                            novelKmers.put(ck, false);
                        }
                    }

                    stretchNum++;
                }
            }
        }

        log.info("Num stretches: {}", stretchNum);
    }
}
