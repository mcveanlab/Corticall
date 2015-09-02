package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.EdgeReversedGraph;
import org.jgrapht.traverse.DepthFirstIterator;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexJunctionsRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksMap;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

public class GenotypeGraphUtils {
    private GenotypeGraphUtils() {}

    private static Set<AnnotatedVertex> predecessorsOf(DirectedGraph<AnnotatedVertex, AnnotatedEdge> sg, String sk) {
        Set<AnnotatedVertex> avs = new HashSet<AnnotatedVertex>();

        DirectedGraph<AnnotatedVertex, AnnotatedEdge> sgr = new EdgeReversedGraph<AnnotatedVertex, AnnotatedEdge>(sg);

        DepthFirstIterator<AnnotatedVertex, AnnotatedEdge> dfs = new DepthFirstIterator<AnnotatedVertex, AnnotatedEdge>(sgr, new AnnotatedVertex(sk));
        while (dfs.hasNext()) {
            avs.add(dfs.next());
        }

        return avs;
    }

    private static Set<AnnotatedVertex> successorsOf(DirectedGraph<AnnotatedVertex, AnnotatedEdge> sg, String sk) {
        Set<AnnotatedVertex> avs = new HashSet<AnnotatedVertex>();

        DepthFirstIterator<AnnotatedVertex, AnnotatedEdge> dfs = new DepthFirstIterator<AnnotatedVertex, AnnotatedEdge>(sg, new AnnotatedVertex(sk));
        while (dfs.hasNext()) {
            avs.add(dfs.next());
        }

        return avs;
    }

    public static DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfsGraph(String stretch, CortexGraph clean, CortexGraph dirty, CortexLinksMap links, boolean beAggressive, Map<CortexKmer, Boolean> novelKmers) {
        CortexGraph GRAPH = clean;
        CortexGraph GRAPH_RAW = dirty;
        boolean AGGRESSIVE = beAggressive;

        // first, explore each color and bring the local subgraphs into memory
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> sg0 = new DefaultDirectedGraph<AnnotatedVertex, AnnotatedEdge>(AnnotatedEdge.class);
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> sg1 = new DefaultDirectedGraph<AnnotatedVertex, AnnotatedEdge>(AnnotatedEdge.class);
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> sg2 = new DefaultDirectedGraph<AnnotatedVertex, AnnotatedEdge>(AnnotatedEdge.class);

        for (int i = 0; i <= stretch.length() - GRAPH.getKmerSize(); i++) {
            String sk = stretch.substring(i, i + GRAPH.getKmerSize());
            CortexKmer ck = new CortexKmer(sk);

            if (novelKmers.containsKey(ck) && !sg0.containsVertex(new AnnotatedVertex(sk, true))) {
                for (boolean goForward : Arrays.asList(true, false)) {
                    TraversalStopper<AnnotatedVertex, AnnotatedEdge> stopper = new AbstractTraversalStopper<AnnotatedVertex, AnnotatedEdge>() {
                        private int goalSize = 0;
                        private int goalDepth = 0;

                        private void reset() {
                            goalSize = 0;
                            goalDepth = 0;
                        }

                        private boolean isLowComplexity(CortexRecord cr) {
                            byte edges[][] = cr.getEdgesAsBytes();

                            int numEdges = 0;
                            for (int e = 0; e < 8; e++) {
                                if (edges[0][e] != '.') {
                                    numEdges++;
                                }
                            }

                            return numEdges > 6;
                        }

                        public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int depth, int size) {
                            /* TODO: this doesn't work if there are multiple branches that need to be truncated.  Fix. */

                            if (goalSize == 0 && (cr.getCoverage(1) > 0 || cr.getCoverage(2) > 0)) {
                                goalSize = size;
                                goalDepth = depth;
                            }

                            if (goalSize > 0 && (cr.getCoverage(0) > 10 && cr.getCoverage(1) == 0 && cr.getCoverage(2) == 0)) {
                                goalSize = size;
                                goalDepth = depth;
                            }

                            if (goalSize > 0 && (size >= goalSize + 50 || isLowComplexity(cr))) {
                                reset();

                                return true;
                            }

                            return false;
                        }

                        @Override
                        public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int depth, int size) {
                            //return goalSize > 0 && depth >= goalDepth + 2;

                            if (goalSize > 0 && depth >= goalDepth + 1) {
                                reset();
                                return true;
                            }

                            return false;
                        }

                        @Override
                        public int maxJunctionsAllowed() {
                            return 10;
                        }
                    };

                    Graphs.addGraph(sg0, CortexUtils.dfs(GRAPH, GRAPH_RAW, sk, 0, null, stopper, goForward));
                }
            }
        }

        //Set<AnnotatedVertex> predecessorList = getPredecessorList(sg0);
        //Set<AnnotatedVertex> successorList = getSuccessorList(sg0);

        //
        Set<AnnotatedVertex> predecessorList = new HashSet<AnnotatedVertex>();
        Set<AnnotatedVertex> successorList = new HashSet<AnnotatedVertex>();

        for (AnnotatedVertex av : sg0.vertexSet()) {
            if (sg0.outDegreeOf(av) != 1) {
                predecessorList.add(av);
            }

            if (sg0.inDegreeOf(av) != 1) {
                successorList.add(av);
            }
        }
        //

        //
        for (int c = 1; c <= 2; c++) {
            DirectedGraph<AnnotatedVertex, AnnotatedEdge> sg = (c == 1) ? sg1 : sg2;

            for (boolean goForward : Arrays.asList(true, false)) {
                Set<AnnotatedVertex> psList = goForward ? predecessorList : successorList;

                for (AnnotatedVertex ak : psList) {
                    TraversalStopper<AnnotatedVertex, AnnotatedEdge> stopper = new AbstractTraversalStopper<AnnotatedVertex, AnnotatedEdge>() {
                        private boolean isLowComplexity(CortexRecord cr) {
                            byte edges[][] = cr.getEdgesAsBytes();

                            int numEdges = 0;
                            for (int e = 0; e < 8; e++) {
                                if (edges[0][e] != '.') {
                                    numEdges++;
                                }
                            }

                            return numEdges > 6;
                        }
                        @Override
                        public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions, int size) {
                            String fw = cr.getKmerAsString();
                            String rc = SequenceUtils.reverseComplement(fw);

                            return size > 1 && (g.containsVertex(new AnnotatedVertex(fw)) || g.containsVertex(new AnnotatedVertex(rc)) || g.containsVertex(new AnnotatedVertex(fw, true)) || g.containsVertex(new AnnotatedVertex(rc, true)));
                        }

                        @Override
                        public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions, int size) {
                            return size > 100 || junctions >= maxJunctionsAllowed() || isLowComplexity(cr);
                        }

                        @Override
                        public int maxJunctionsAllowed() {
                            return 3;
                        }
                    };

                    DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs = CortexUtils.dfs(GRAPH, GRAPH_RAW, ak.getKmer(), c, sg0, stopper, goForward);
                    if (stopper.anyTraversalSucceeded()) {
                        Graphs.addGraph(sg, dfs);
                    }
                }
            }
        }
        //

        // Now, combine them all into an annotated graph
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> ag = new DefaultDirectedGraph<AnnotatedVertex, AnnotatedEdge>(AnnotatedEdge.class);
        addGraph(ag, sg0, 0, novelKmers);
        //addGraph(ag, sg1, 1, novelKmers);
        //addGraph(ag, sg2, 2, novelKmers);

        for (AnnotatedVertex ava : ag.vertexSet()) {
            if (predecessorList.contains(ava)) { ava.setFlag("predecessor"); }
            if (successorList.contains(ava)) { ava.setFlag("successor"); }
        }

        return ag;
    }

    private static Set<AnnotatedVertex> getPredecessorList(DirectedGraph<AnnotatedVertex, AnnotatedEdge> sg0) {
        Set<AnnotatedVertex> predecessorList = new HashSet<AnnotatedVertex>();

        for (AnnotatedVertex av : sg0.vertexSet()) {
            Set<AnnotatedEdge> oes = sg0.outgoingEdgesOf(av);

            boolean hasOutgoingNovelEdge = false;
            Set<AnnotatedEdge> outgoingParentalEdges = new HashSet<AnnotatedEdge>();

            for (AnnotatedEdge oe : oes) {
                if (oe.isPresent(0) && oe.isAbsent(1) && oe.isAbsent(2)) {
                    hasOutgoingNovelEdge = true;
                } else if (oe.isPresent(1) || oe.isPresent(2)) {
                    outgoingParentalEdges.add(oe);
                }
            }

            if (hasOutgoingNovelEdge && outgoingParentalEdges.size() > 0) {
                for (AnnotatedEdge oe : outgoingParentalEdges) {
                    AnnotatedVertex v = sg0.getEdgeTarget(oe);

                    DepthFirstIterator<AnnotatedVertex, AnnotatedEdge> dfs = new DepthFirstIterator<AnnotatedVertex, AnnotatedEdge>(sg0, v);
                    while (dfs.hasNext()) {
                        AnnotatedVertex nv = dfs.next();

                        if (sg0.outDegreeOf(nv) == 0) {
                            predecessorList.add(nv);
                        }
                    }
                }
            }
        }

        return predecessorList;
    }

    private static Set<AnnotatedVertex> getSuccessorList(DirectedGraph<AnnotatedVertex, AnnotatedEdge> sg0) {
        Set<AnnotatedVertex> successorList = new HashSet<AnnotatedVertex>();

        for (AnnotatedVertex av : sg0.vertexSet()) {
            Set<AnnotatedEdge> ies = sg0.incomingEdgesOf(av);

            boolean hasIncomingNovelEdge = false;
            Set<AnnotatedEdge> incomingParentalEdges = new HashSet<AnnotatedEdge>();

            for (AnnotatedEdge ie : ies) {
                if (ie.isPresent(0) && ie.isAbsent(1) && ie.isAbsent(2)) {
                    hasIncomingNovelEdge = true;
                } else if (ie.isPresent(1) || ie.isPresent(2)) {
                    incomingParentalEdges.add(ie);
                }
            }

            if (hasIncomingNovelEdge && incomingParentalEdges.size() > 0) {
                for (AnnotatedEdge ie : incomingParentalEdges) {
                    AnnotatedVertex v = sg0.getEdgeSource(ie);

                    DepthFirstIterator<AnnotatedVertex, AnnotatedEdge> rdfs = new DepthFirstIterator<AnnotatedVertex, AnnotatedEdge>(new EdgeReversedGraph<AnnotatedVertex, AnnotatedEdge>(sg0), v);
                    while (rdfs.hasNext()) {
                        AnnotatedVertex pv = rdfs.next();

                        if (sg0.inDegreeOf(pv) == 0) {
                            successorList.add(pv);
                        }
                    }
                }
            }
        }

        return successorList;
    }

    private static String formatAttributes(Map<String, Object> attrs) {
        List<String> attrArray = new ArrayList<String>();

        for (String attr : attrs.keySet()) {
            String value = attrs.get(attr).toString();

            attrArray.add(attr + "=\"" + value + "\"");
        }

        return Joiner.on(" ").join(attrArray);
    }

    public static void printGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, File out, boolean withText, boolean withPdf) {
        try {
            File f = new File(out.getAbsolutePath() + ".dot");
            File p = new File(out.getAbsolutePath() + ".pdf");

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

                if (!withText || g.inDegreeOf(v) == 0 || g.outDegreeOf(v) == 0) {
                    attrs.put("label", "");
                }

                if (v.isNovel()) {
                    attrs.put("color", "red");
                    attrs.put("fillcolor", "red");
                    attrs.put("shape", withText ? "box" : "circle");
                }

                if (v.flagIsSet("start") || v.flagIsSet("end")) {
                    attrs.put("label", v.getKmer());
                    attrs.put("fillcolor", v.flagIsSet("start") ? "orange" : "purple");
                    attrs.put("shape", "rect");
                }

                o.println(indent + "\"" + v.getKmer() + "\" [ " + formatAttributes(attrs) + " ];");
            }

            String[] colors = new String[]{"red", "blue", "green"};

            for (AnnotatedEdge e : g.edgeSet()) {
                String s = g.getEdgeSource(e).getKmer();
                String t = g.getEdgeTarget(e).getKmer();

                for (int c = 0; c < 3; c++) {
                    if (e.isPresent(c)) {
                        Map<String, Object> attrs = new TreeMap<String, Object>();
                        attrs.put("color", colors[c]);
                        attrs.put("weight", g.getEdgeWeight(e));
                        attrs.put("label", g.getEdgeWeight(e));

                        o.println(indent + "\"" + s + "\" -> \"" + t + "\" [ " + formatAttributes(attrs) + " ];");
                    }
                }
            }

            o.println("}");

            o.close();

            if (withPdf) {
                Runtime.getRuntime().exec("dot -Tpdf -o" + p.getAbsolutePath() + " " + f.getAbsolutePath());
            }
        } catch (FileNotFoundException e) {
            throw new IndianaException("File not found", e);
        } catch (IOException e) {
            throw new IndianaException("IO exception", e);
        }
    }

    public static void addGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int color, Map<CortexKmer, Boolean> novelKmers) {
        for (AnnotatedVertex v : g.vertexSet()) {
            AnnotatedVertex av = new AnnotatedVertex(v.getKmer(), novelKmers.containsKey(new CortexKmer(v.getKmer())));

            a.addVertex(av);
        }

        for (AnnotatedEdge e : g.edgeSet()) {
            AnnotatedVertex s0 = g.getEdgeSource(e);
            AnnotatedVertex s1 = g.getEdgeTarget(e);

            CortexKmer ck0 = new CortexKmer(s0.getKmer());
            CortexKmer ck1 = new CortexKmer(s1.getKmer());

            AnnotatedVertex a0 = new AnnotatedVertex(s0.getKmer(), novelKmers.containsKey(ck0));
            AnnotatedVertex a1 = new AnnotatedVertex(s1.getKmer(), novelKmers.containsKey(ck1));

            if (!a.containsEdge(a0, a1)) {
                a.addEdge(a0, a1, new AnnotatedEdge(e.getPresence()));
            }

            //a.getEdge(a0, a1).set(color, true);
        }
    }

    public static DirectedGraph<AnnotatedVertex, AnnotatedEdge> copyGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a) {
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

    public static DirectedGraph<AnnotatedVertex, AnnotatedEdge> simplifyGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, boolean combine) {
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

                while (nextVertex != null && b.outDegreeOf(thisVertex) == 1 && b.inDegreeOf(nextVertex) == 1 && (combine || thisVertex.isNovel() == nextVertex.isNovel())) {
                    AnnotatedVertex sv = new AnnotatedVertex(thisVertex.getKmer() + nextVertex.getKmer().charAt(nextVertex.getKmer().length() - 1), thisVertex.isNovel());

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

    public static Map<String, VariantInfo> loadVariantInfoMap(File bedFile) {
        Map<String, VariantInfo> allVariants = new HashMap<String, VariantInfo>();

        if (bedFile != null) {
            TableReader bed = new TableReader(bedFile, "vchr", "vstart", "vstop", "info");

            for (Map<String, String> be : bed) {
                VariantInfo vi = new VariantInfo();

                vi.vchr = be.get("vchr");
                vi.vstart = Integer.valueOf(be.get("vstart"));
                vi.vstop = Integer.valueOf(be.get("vstop"));

                String[] kvpairs = be.get("info").split(";");
                for (String kvpair : kvpairs) {
                    String[] kv = kvpair.split("=");

                    if (kv[0].equals("id")) {
                        vi.variantId = kv[1];
                    }
                    if (kv[0].equals("type")) {
                        vi.type = kv[1];
                    }
                    if (kv[0].equals("denovo")) {
                        vi.denovo = kv[1];
                    }
                    if (kv[0].equals("nahr")) {
                        vi.nahr = kv[1];
                    }
                    if (kv[0].equals("gcindex")) {
                        vi.gcindex = Integer.valueOf(kv[1]);
                    }
                    if (kv[0].equals("ref")) {
                        vi.ref = kv[1];
                    }
                    if (kv[0].equals("alt")) {
                        vi.alt = kv[1];
                    }
                    if (kv[0].equals("left")) {
                        vi.leftFlank = kv[1];
                    }
                    if (kv[0].equals("right")) {
                        vi.rightFlank = kv[1];
                    }
                }

                allVariants.put(vi.variantId, vi);
            }
        }
        return allVariants;
    }

    public static Map<CortexKmer, VariantInfo> loadNovelKmerMap(File novelKmerMap, File bedFile) {
        Map<String, VariantInfo> allVariants = loadVariantInfoMap(bedFile);

        Map<CortexKmer, VariantInfo> vis = new HashMap<CortexKmer, VariantInfo>();

        if (novelKmerMap != null) {
            TableReader tr = new TableReader(novelKmerMap);

            for (Map<String, String> te : tr) {
                VariantInfo vi = allVariants.get(te.get("variantId"));
                CortexKmer kmer = new CortexKmer(te.get("kmer"));

                vis.put(kmer, vi);
            }
        }

        return vis;
    }
}
