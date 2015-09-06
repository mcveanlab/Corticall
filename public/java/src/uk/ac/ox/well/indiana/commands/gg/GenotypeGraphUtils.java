package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
import htsjdk.samtools.util.Interval;
import org.jgrapht.DirectedGraph;
import org.jgrapht.GraphPath;
import org.jgrapht.Graphs;
import org.jgrapht.alg.DijkstraShortestPath;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.EdgeReversedGraph;
import org.jgrapht.traverse.DepthFirstIterator;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.math.MoreMathUtils;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

public class GenotypeGraphUtils {
    private GenotypeGraphUtils() {}

    private static int numRemainingLocalNovelKmers(Map<String, Boolean> localNovelKmers) {
        int numRemainingLocalNovelKmers = 0;

        for (String sk : localNovelKmers.keySet()) {
            if (localNovelKmers.get(sk)) {
                numRemainingLocalNovelKmers++;
            }
        }

        return numRemainingLocalNovelKmers;
    }

    public static DirectedGraph<AnnotatedVertex, AnnotatedEdge> loadLocalSubgraph(String stretch, CortexGraph clean, CortexGraph dirty, Map<CortexKmer, Boolean> novelKmers) {
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> sg0 = new DefaultDirectedGraph<AnnotatedVertex, AnnotatedEdge>(AnnotatedEdge.class);
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> sg1 = new DefaultDirectedGraph<AnnotatedVertex, AnnotatedEdge>(AnnotatedEdge.class);
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> sg2 = new DefaultDirectedGraph<AnnotatedVertex, AnnotatedEdge>(AnnotatedEdge.class);

        Map<String, Boolean> localNovelKmers = new HashMap<String, Boolean>();
        for (int i = 0; i <= stretch.length() - clean.getKmerSize(); i++) {
            String sk = stretch.substring(i, i + clean.getKmerSize());
            CortexKmer ck = new CortexKmer(sk);

            if (novelKmers.containsKey(ck)) {
                localNovelKmers.put(sk, true);
            }
        }

        do {
            for (String sk : localNovelKmers.keySet()) {
                if (localNovelKmers.get(sk)) {
                    DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs = CortexUtils.dfs(clean, dirty, sk, 0, null, ChildTraversalStopper.class);

                    if (dfs != null) {
                        Graphs.addGraph(sg0, dfs);
                    }

                    localNovelKmers.put(sk, false);
                }
            }

            Set<AnnotatedVertex> predecessorList = new HashSet<AnnotatedVertex>();
            Set<AnnotatedVertex> successorList = new HashSet<AnnotatedVertex>();

            for (AnnotatedVertex av : sg0.vertexSet()) {
                if (!novelKmers.containsKey(new CortexKmer(av.getKmer()))) {
                    if (sg0.outDegreeOf(av) > 1) {
                        Set<AnnotatedEdge> novelEdges = new HashSet<AnnotatedEdge>();
                        Set<AnnotatedEdge> parentalEdges = new HashSet<AnnotatedEdge>();

                        for (AnnotatedEdge ae : sg0.outgoingEdgesOf(av)) {
                            if (ae.isPresent(0) && ae.isAbsent(1) && ae.isAbsent(2)) {
                                novelEdges.add(ae);
                            } else if (ae.isPresent(1) || ae.isPresent(2)) {
                                parentalEdges.add(ae);
                            }
                        }

                        if (novelEdges.size() > 0 && parentalEdges.size() > 0) {
                            for (AnnotatedEdge pe : parentalEdges) {
                                AnnotatedVertex tv = sg0.getEdgeTarget(pe);

                                DepthFirstIterator<AnnotatedVertex, AnnotatedEdge> dfi = new DepthFirstIterator<AnnotatedVertex, AnnotatedEdge>(sg0, tv);

                                while (dfi.hasNext()) {
                                    AnnotatedVertex nv = dfi.next();

                                    if (sg0.outDegreeOf(nv) == 0) {
                                        predecessorList.add(nv);
                                    }
                                }
                            }

                        }
                    }

                    if (sg0.inDegreeOf(av) > 1) {
                        Set<AnnotatedEdge> novelEdges = new HashSet<AnnotatedEdge>();
                        Set<AnnotatedEdge> parentalEdges = new HashSet<AnnotatedEdge>();

                        for (AnnotatedEdge ae : sg0.incomingEdgesOf(av)) {
                            if (ae.isPresent(0) && ae.isAbsent(1) && ae.isAbsent(2)) {
                                novelEdges.add(ae);
                            } else if (ae.isPresent(1) || ae.isPresent(2)) {
                                parentalEdges.add(ae);
                            }
                        }

                        if (novelEdges.size() > 0 && parentalEdges.size() > 0) {
                            for (AnnotatedEdge pe : parentalEdges) {
                                AnnotatedVertex sv = sg0.getEdgeSource(pe);

                                DepthFirstIterator<AnnotatedVertex, AnnotatedEdge> dfi = new DepthFirstIterator<AnnotatedVertex, AnnotatedEdge>(new EdgeReversedGraph<AnnotatedVertex, AnnotatedEdge>(sg0), sv);

                                while (dfi.hasNext()) {
                                    AnnotatedVertex pv = dfi.next();

                                    if (sg0.inDegreeOf(pv) == 0) {
                                        successorList.add(pv);
                                    }
                                }
                            }

                        }
                    }
                }
            }

            for (AnnotatedVertex ava : sg0.vertexSet()) {
                if (predecessorList.contains(ava)) {
                    ava.setFlag("predecessor");
                }
                if (successorList.contains(ava)) {
                    ava.setFlag("successor");
                }
            }

            //
            for (int c = 1; c <= 2; c++) {
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> sg = (c == 1) ? sg1 : sg2;

                for (boolean goForward : Arrays.asList(true, false)) {
                    Set<AnnotatedVertex> psList = goForward ? predecessorList : successorList;

                    for (AnnotatedVertex ak : psList) {
                        DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs = CortexUtils.dfs(clean, dirty, ak.getKmer(), c, sg0, ParentTraversalStopper.class, 0, goForward);

                        if (dfs != null) {
                            Graphs.addGraph(sg, dfs);
                        }
                    }
                }

                for (AnnotatedVertex av : sg.vertexSet()) {
                    CortexKmer ck = new CortexKmer(av.getKmer());

                    if (novelKmers.containsKey(ck) && !localNovelKmers.containsKey(av.getKmer())) {
                        localNovelKmers.put(av.getKmer(), true);
                    }
                }
            }
            //

        } while (numRemainingLocalNovelKmers(localNovelKmers) > 0);

        // Now, combine them all into an annotated graph
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> ag = new DefaultDirectedGraph<AnnotatedVertex, AnnotatedEdge>(AnnotatedEdge.class);
        addGraph(ag, sg0, 0, novelKmers);
        addGraph(ag, sg1, 1, novelKmers);
        addGraph(ag, sg2, 2, novelKmers);

        /*
        for (AnnotatedVertex ava : ag.vertexSet()) {
            if (predecessorList.contains(ava)) {
                ava.setFlag("predecessor");
            }
            if (successorList.contains(ava)) {
                ava.setFlag("successor");
            }
        }
        */

        return ag;
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

    private static DirectedGraph<AnnotatedVertex, AnnotatedEdge> removeOtherColors(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, int colorToRetain) {
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

    private static List<String> getNovelStretch(String stretch, Map<CortexKmer, Boolean> novelKmers) {
        boolean[] novel = new boolean[stretch.length()];
        int kmerLength = novelKmers.keySet().iterator().next().length();

        for (int i = 0; i <= stretch.length() - kmerLength; i++) {
            String sk = stretch.substring(i, i + kmerLength);
            CortexKmer ck = new CortexKmer(sk);

            novel[i] = novelKmers.containsKey(ck);
        }

        List<String> novelStretches = new ArrayList<String>();

        StringBuilder novelStretchBuilder = new StringBuilder();

        for (int i = 0; i <= stretch.length() - kmerLength; i++) {
            String sk = stretch.substring(i, i + kmerLength);
            boolean isNovel = novel[i];

            if (!isNovel) {
                if (novelStretchBuilder.length() > 0) {
                    novelStretches.add(novelStretchBuilder.toString());
                    novelStretchBuilder = new StringBuilder();
                }
            } else {
                if (novelStretchBuilder.length() == 0) {
                    novelStretchBuilder.append(sk);
                } else {
                    novelStretchBuilder.append(sk.charAt(sk.length() - 1));
                }
            }
        }

        if (novelStretchBuilder.length() > 0) {
            novelStretches.add(novelStretchBuilder.toString());
        }

        return novelStretches;
    }

    private static String linearizePath(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, GraphPath<AnnotatedVertex, AnnotatedEdge> p) {
        StringBuilder sb = new StringBuilder(p.getStartVertex().getKmer());

        for (AnnotatedEdge ae : p.getEdgeList()) {
            AnnotatedVertex av = a.getEdgeTarget(ae);

            sb.append(av.getKmer().charAt(av.getKmer().length() - 1));
        }

        return sb.toString();
    }

    public static PathInfo computeBestMinWeightPath(CortexGraph clean, CortexGraph dirty, DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, int color, String stretch, Map<CortexKmer, Boolean> novelKmers) {
        class AnnotateStartsAndEnds {
            private int color;
            private String stretch;
            private Map<CortexKmer, Boolean> novelKmers;
            private DirectedGraph<AnnotatedVertex, AnnotatedEdge> b;
            private Set<AnnotatedVertex> candidateEnds = new HashSet<AnnotatedVertex>();
            private Set<AnnotatedVertex> candidateStarts = new HashSet<AnnotatedVertex>();

            public AnnotateStartsAndEnds(int color, String stretch, Map<CortexKmer, Boolean> novelKmers, DirectedGraph<AnnotatedVertex, AnnotatedEdge> b) {
                this.color = color;
                this.stretch = stretch;
                this.novelKmers = novelKmers;
                this.b = b;
            }

            public Set<AnnotatedVertex> getCandidateEnds() {
                return candidateEnds;
            }

            public Set<AnnotatedVertex> getCandidateStarts() { return candidateStarts; }

            public AnnotateStartsAndEnds invoke() {
                for (AnnotatedVertex av : b.vertexSet()) {
                    if (b.outDegreeOf(av) > 1) {
                        Set<AnnotatedEdge> aes = b.outgoingEdgesOf(av);

                        Set<AnnotatedEdge> childEdges = new HashSet<AnnotatedEdge>();
                        Set<AnnotatedEdge> parentEdges = new HashSet<AnnotatedEdge>();

                        for (AnnotatedEdge ae : aes) {
                            if (ae.isPresent(0) && ae.isAbsent(1) && ae.isAbsent(2)) {
                                childEdges.add(ae);
                            }

                            if (ae.isPresent(color)) {
                                parentEdges.add(ae);
                            }
                        }

                        if (childEdges.size() > 0 && parentEdges.size() > 0) {
                            av.setFlag("start");

                            candidateStarts.add(av);
                        }
                    }

                    if (b.inDegreeOf(av) > 1) {
                        Set<AnnotatedEdge> aes = b.incomingEdgesOf(av);

                        Set<AnnotatedEdge> childEdges = new HashSet<AnnotatedEdge>();
                        Set<AnnotatedEdge> parentEdges = new HashSet<AnnotatedEdge>();

                        for (AnnotatedEdge ae : aes) {
                            if (ae.isPresent(0) && ae.isAbsent(1) && ae.isAbsent(2)) {
                                childEdges.add(ae);
                            }

                            if (ae.isPresent(color)) {
                                parentEdges.add(ae);
                            }
                        }

                        if (childEdges.size() > 0 && parentEdges.size() > 0) {
                            av.setFlag("end");

                            candidateEnds.add(av);
                        }
                    }
                }

                return this;
            }
        }

        DirectedGraph<AnnotatedVertex, AnnotatedEdge> b = copyGraph(a);

        AnnotateStartsAndEnds annotateStartsAndEnds = new AnnotateStartsAndEnds(color, stretch, novelKmers, b).invoke();
        Set<AnnotatedVertex> candidateStarts = annotateStartsAndEnds.getCandidateStarts();
        Set<AnnotatedVertex> candidateEnds = annotateStartsAndEnds.getCandidateEnds();

        DirectedGraph<AnnotatedVertex, AnnotatedEdge> b0 = removeOtherColors(b, 0);
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> bc = removeOtherColors(b, color);

        List<String> novelStretches = getNovelStretch(stretch, novelKmers);

        double minPl0 = Double.MAX_VALUE, minPlc = Double.MAX_VALUE;
        String minLp0 = "", minLpc = "";
        String start = "", stop = "";

        Set<AnnotatedVertex> sought = new HashSet<AnnotatedVertex>();
        sought.add(new AnnotatedVertex("GATAATATATTTTGTAGAATATTTTTATTATTTTTAATTAACATGAA"));
        sought.add(new AnnotatedVertex("AAGAATTTTTTTACACTTTTATTTTCGTTATTTGTATCTTCATTATT"));
        sought.add(new AnnotatedVertex("TTATATGTGGAAGAATTATGTATAATATTATTATTATCATTATTATT"));
        sought.add(new AnnotatedVertex("TTTTTTTGCATATTTAATTGATCCATTTTATTATAATCCTTATCATT"));

        for (AnnotatedVertex sv : candidateStarts) {
            for (AnnotatedVertex ev : candidateEnds) {
                String lp0, lpc;

                if (b0.containsVertex(sv) && b0.containsVertex(ev) && bc.containsVertex(sv) && bc.containsVertex(ev)) {
                    do {
                        DijkstraShortestPath<AnnotatedVertex, AnnotatedEdge> dsp0 = new DijkstraShortestPath<AnnotatedVertex, AnnotatedEdge>(b0, sv, ev);
                        DijkstraShortestPath<AnnotatedVertex, AnnotatedEdge> dspc = new DijkstraShortestPath<AnnotatedVertex, AnnotatedEdge>(bc, sv, ev);

                        GraphPath<AnnotatedVertex, AnnotatedEdge> p0 = dsp0.getPath();
                        GraphPath<AnnotatedVertex, AnnotatedEdge> pc = dspc.getPath();

                        lp0 = p0 == null ? "" : linearizePath(b0, p0);
                        lpc = pc == null ? "" : linearizePath(bc, pc);

                        if (lp0.contains(sv.getKmer()) && lp0.contains(ev.getKmer()) && lpc.contains(sv.getKmer()) && lpc.contains(ev.getKmer())) {
                            boolean stretchesArePresent = true;
                            for (String novelStretch : novelStretches) {
                                if (!lp0.contains(novelStretch)) {
                                    stretchesArePresent = false;
                                    break;
                                }
                            }

                            if (novelStretches.size() > 0 && stretchesArePresent) {
                                double pl0 = dsp0.getPathLength();
                                double plc = dspc.getPathLength();

                                if (pl0 < minPl0 && plc < minPlc) {
                                    minPl0 = pl0;
                                    minLp0 = lp0;

                                    minPlc = plc;
                                    minLpc = lpc;

                                    start = sv.getKmer();
                                    stop = ev.getKmer();
                                }
                            }
                        }

                        if (lp0.equals(lpc) && !lp0.isEmpty() && lp0.length() > clean.getKmerSize() + 1) {
                            for (int i = 1; i < lp0.length() - clean.getKmerSize(); i++) {
                                AnnotatedVertex as = new AnnotatedVertex(lp0.substring(i, i + clean.getKmerSize()));

                                b0.removeVertex(as);
                            }
                        }
                    } while (lp0.equals(lpc) && !lp0.isEmpty() && lp0.length() > clean.getKmerSize() + 1);
                }
            }
        }

        return new PathInfo(start, stop, minLp0, minLpc);
    }

    public static GraphicalVariantContext callVariant(CortexGraph clean, CortexGraph dirty, PathInfo p, int color, String stretch, Map<CortexKmer, Boolean> novelKmers, KmerLookup kl) {
        // Compute paths
        //PathInfo p = computeBestMinWeightPath(a, color, stretch, novelKmers);

        // Trim back to reference and variant alleles
        int s, e0 = p.child.length() - 1, e1 = p.parent.length() - 1;

        for (s = 0; s < (p.child.length() < p.parent.length() ? p.child.length() : p.parent.length()) && p.child.charAt(s) == p.parent.charAt(s); s++) {
        }

        while (e0 > s && e1 > s && p.child.charAt(e0) == p.parent.charAt(e1)) {
            e0--;
            e1--;
        }

        String parentalAllele = p.parent == null || p.parent.equals("") || p.child.equals(p.parent) ? "A" : p.parent.substring(s, e1 + 1);
        String childAllele = p.child == null || p.child.equals("") || p.child.equals(p.parent) ? "N" : p.child.substring(s, e0 + 1);

        int e = s + parentalAllele.length() - 1;

        // Decide if the event is actually a GC or NAHR event
        boolean hasRecombs = hasRecombinations(clean, dirty, stretch);
        boolean isChimeric = isChimeric(stretch, kl);

        List<Set<Interval>> alignment = kl.align(CortexUtils.getSeededStretchLeft(clean, p.start, color, false) + p.parent + CortexUtils.getSeededStretchRight(clean, p.stop, color, false));
        List<Set<Interval>> anovel = kl.align(stretch);

        // Build the GVC
        GraphicalVariantContext gvc = new GraphicalVariantContext()
                .attribute(color, "start", s)
                .attribute(color, "stop", e)
                .attribute(color, "parentalAllele", parentalAllele)
                .attribute(color, "childAllele", childAllele)
                .attribute(color, "parentalStretch", p.parent)
                .attribute(color, "childStretch", p.child)
                .attribute(color, "event", "unknown")
                .attribute(color, "traversalStatus", "complete")
                .attribute(color, "parentalPathAlignment", alignment)
                .attribute(color, "novelStretchAlignment", anovel)
                .attribute(color, "haplotypeBackground", color);

        if (childAllele.equals("N")) {
            if (hasRecombs) {
                gvc.attribute(color, "event", "GC");
            } else if (isChimeric) {
                gvc.attribute(color, "event", "NAHR");
            } else {
                gvc.attribute(color, "traversalStatus", "incomplete");
            }
        } else {
            if (parentalAllele.length() == childAllele.length()) {
                if (parentalAllele.length() == 1) {
                    gvc.attribute(color, "event", "SNP");
                } else if (SequenceUtils.reverseComplement(parentalAllele).equals(childAllele)) {
                    gvc.attribute(color, "event", "INV");
                }
            } else if (parentalAllele.length() == 1 && childAllele.length() > 1) {
                gvc.attribute(color, "event", "INS");

                String childAlleleTrimmed = childAllele;
                if (parentalAllele.length() == 1 && childAllele.charAt(0) == parentalAllele.charAt(0)) {
                    childAlleleTrimmed = childAllele.substring(1, childAllele.length());
                } else if (parentalAllele.length() == 1 && childAllele.charAt(childAllele.length() - 1) == parentalAllele.charAt(0)) {
                    childAlleleTrimmed = childAllele.substring(0, childAllele.length() - 1);
                }

                String repUnit = "";
                boolean isStrExp = false;
                for (int repLength = 2; repLength <= 5 && repLength <= childAlleleTrimmed.length(); repLength++) {
                    repUnit = childAlleleTrimmed.substring(0, repLength);

                    boolean fits = true;
                    for (int i = 0; i < childAlleleTrimmed.length(); i += repLength) {
                        fits &= (i + repLength <= childAllele.length()) && childAllele.substring(i, i + repLength).equals(repUnit);
                    }

                    if (fits &&
                            (p.parent.substring(s - repLength, s).equals(repUnit) ||
                                    p.parent.substring(s, s + repLength).equals(repUnit))) {
                        gvc.attribute(color, "event", "STR_EXP");
                        gvc.attribute(color, "repeatingUnit", repUnit);
                        isStrExp = true;
                    }
                }

                if (!isStrExp && childAlleleTrimmed.length() >= 10 && childAlleleTrimmed.length() <= 50 &&
                        ((s - childAlleleTrimmed.length() >= 0 && p.parent.substring(s - childAlleleTrimmed.length(), s).equals(childAlleleTrimmed)) ||
                                (s + childAlleleTrimmed.length() <= p.parent.length() && p.parent.substring(s, s + childAlleleTrimmed.length()).equals(childAlleleTrimmed)))) {
                    gvc.attribute(color, "event", "TD");
                    gvc.attribute(color, "repeatingUnit", repUnit);
                }
            } else if (parentalAllele.length() > 1 && childAllele.length() == 1) {
                gvc.attribute(color, "event", "DEL");

                String parentalAlleleTrimmed = parentalAllele;
                if (childAllele.length() == 1 && parentalAllele.charAt(0) == childAllele.charAt(0)) {
                    parentalAlleleTrimmed = parentalAllele.substring(1, parentalAllele.length());
                } else if (childAllele.length() == 1 && parentalAllele.charAt(parentalAllele.length() - 1) == childAllele.charAt(0)) {
                    parentalAlleleTrimmed = parentalAllele.substring(0, parentalAllele.length() - 1);
                }

                String repUnit = "";
                for (int repLength = 2; repLength <= 5 && repLength <= parentalAlleleTrimmed.length(); repLength++) {
                    repUnit = parentalAlleleTrimmed.substring(0, repLength);

                    boolean fits = true;
                    for (int i = 0; i < parentalAlleleTrimmed.length(); i += repLength) {
                        fits &= (i + repLength <= parentalAllele.length()) && parentalAllele.substring(i, i + repLength).equals(repUnit);
                    }

                    if (fits &&
                            (p.child.substring(s - repLength, s).equals(repUnit) ||
                                    p.child.substring(s + parentalAlleleTrimmed.length(), s + parentalAlleleTrimmed.length() + repLength).equals(repUnit))) {
                        gvc.attribute(color, "event", "STR_CON");
                        gvc.attribute(color, "repeatingUnit", repUnit);
                    }
                }
            }
        }

        return gvc;
    }

    private static boolean hasRecombinations(CortexGraph clean, CortexGraph dirty, String stretch) {
        StringBuilder inherit = new StringBuilder();
        for (int i = 0; i <= stretch.length() - clean.getKmerSize(); i++) {
            String sk = stretch.substring(i, i + clean.getKmerSize());
            CortexKmer ck = new CortexKmer(sk);

            CortexRecord cr = clean.findRecord(ck);

            if (cr != null) {
                int cov0 = cr.getCoverage(0);
                int cov1 = cr.getCoverage(1);
                int cov2 = cr.getCoverage(2);

                if (cov0 > 0 && cov1 == 0 && cov2 == 0) {
                    inherit.append("C");
                } else if (cov0 > 0 && cov1 > 0 && cov2 == 0) {
                    inherit.append("M");
                } else if (cov0 > 0 && cov1 == 0 && cov2 > 0) {
                    inherit.append("D");
                } else if (cov0 > 0 && cov1 > 0 && cov2 > 0) {
                    inherit.append("B");
                } else {
                    inherit.append("?");
                }
            }
        }

        for (int i = 0; i < inherit.length(); i++) {
            if (inherit.charAt(i) == 'B') {
                char prevContext = '?';
                char nextContext = '?';

                for (int j = i - 1; j >= 0; j--) {
                    if (inherit.charAt(j) == 'M' || inherit.charAt(j) == 'D') {
                        prevContext = inherit.charAt(j);
                        break;
                    }
                }

                for (int j = i + 1; j < inherit.length(); j++) {
                    if (inherit.charAt(j) == 'M' || inherit.charAt(j) == 'D') {
                        nextContext = inherit.charAt(j);
                        break;
                    }
                }

                char context = '?';
                if (prevContext == nextContext && prevContext != '?') {
                    context = prevContext;
                } else if (prevContext != nextContext) {
                    if (prevContext != '?') {
                        context = prevContext;
                    } else {
                        context = nextContext;
                    }
                }

                inherit.setCharAt(i, context);
            }
        }

        String inheritStr = inherit.toString();
        return inheritStr.matches(".*D+.C+.M+.*") || inheritStr.matches(".*M+.C+.D+.*");
    }

    private static boolean isChimeric(String stretch, KmerLookup kl) {
        List<Set<Interval>> ils = kl.findKmers(stretch);

        StringBuilder sb = new StringBuilder();

        int nextAvailableCode = 0;
        Map<String, Integer> chrCodes = new HashMap<String, Integer>();
        Map<String, Integer> chrCount = new HashMap<String, Integer>();

        for (int i = 0; i < ils.size(); i++) {
            Set<Interval> il = ils.get(i);

            if (il.size() == 1) {
                Interval interval = il.iterator().next();

                ContainerUtils.increment(chrCount, interval.getSequence());

                if (!chrCodes.containsKey(interval.getSequence())) {
                    chrCodes.put(interval.getSequence(), nextAvailableCode);
                    nextAvailableCode++;
                }

                int code = chrCodes.get(interval.getSequence());
                sb.append(code);
            } else {
                sb.append(".");
            }
        }

        return chrCount.size() > 1;
    }

    public static void chooseVariant(GraphicalVariantContext gvc) {
        int[] scores = new int[2];

        for (int c = 1; c <= 2; c++) {
            if (gvc.getAttribute(c, "traversalStatus").equals("complete")) { scores[c-1]++; }
            if (!gvc.getAttribute(c, "event").equals("unknown")) { scores[c-1]++; }
            if (gvc.getAttribute(c, "event").equals("GC") || gvc.getAttribute(c, "event").equals("NAHR")) { scores[c-1]++; }

            gvc.attribute(c, "score", scores[c-1]);
        }

        int bc = MoreMathUtils.whichMax(scores) + 1;

        gvc.addAttributes(0, gvc.getAttributes(bc));
        gvc.attribute(0, "haplotypeBackground", scores[0] == scores[1] ? 0 : bc);
    }

    public static String recordToString(String sk, CortexRecord cr) {
        String kmer = "";
        String cov = "";
        String ed = "";

        if (cr != null) {
            kmer = cr.getKmerAsString();

            boolean fw = sk.equals(kmer);

            if (!fw) {
                kmer = SequenceUtils.reverseComplement(kmer);
            }

            for (int coverage : cr.getCoverages()) {
                cov += " " + coverage;
            }

            for (String edge : cr.getEdgeAsStrings()) {
                ed += " " + (fw ? edge : SequenceUtils.reverseComplement(edge));
            }
        }

        return kmer + " " + cov + " " + ed;
    }

}
