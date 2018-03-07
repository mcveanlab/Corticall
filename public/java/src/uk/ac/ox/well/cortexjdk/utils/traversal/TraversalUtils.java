package uk.ac.ox.well.cortexjdk.utils.traversal;

import org.jgrapht.GraphPath;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.DestinationStopper;

import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.REVERSE;

public class TraversalUtils {
    private TraversalUtils() {}

    /*
    public static Map<Integer, List<List<List<CortexVertex>>>> fill(List<CortexVertex> w, CortexGraph graph, List<CortexLinks> links, Set<Integer> sampleColors, List<Set<Integer>> backgroundColors) {
        Map<Integer, List<List<List<CortexVertex>>>> mss = new HashMap<>();

        for (Set<Integer> cs : backgroundColors) {
            List<List<CortexVertex>> ms = new ArrayList<>();
            List<CortexVertex> m = new ArrayList<>();

            for (int i = 0; i < w.size(); i++) {
                CortexVertex v0 = w.get(i);
                CortexVertex v1 = i + 1 < w.size() ? w.get(i + 1) : null;

                Set<String> remainingNext = convertToStrings(TraversalUtils.getAllNextKmers(v0.getCortexRecord(), !v0.getKmerAsString().equals(v0.getCanonicalKmer().getKmerAsString())).get(c));

                for (int c : cs) {
                    if (v0.getCortexRecord().getCoverage(c) > 0) {
                        m.add(v0);
                        break;
                    }
                }

                if (v1 != null && !remainingNext.contains(v1.getKmerAsString())) {
                    if (m.size() > 0) {
                        ms.add(m);
                    }

                    m = new ArrayList<>();
                }
            }

            if (m.size() > 0) {
                ms.add(m);
            }

            TraversalEngine ef = new TraversalEngineFactory()
                    .traversalColors(cs)
                    //.recruitmentColors(sampleColor)
                    .traversalDirection(FORWARD)
                    .combinationOperator(OR)
                    .stoppingRule(DestinationStopper.class)
                    .maxBranchLength(1000)
                    .graph(graph)
                    .links(links)
                    .make();

            TraversalEngine er = new TraversalEngineFactory()
                    .traversalColors(cs)
                    //.recruitmentColors(sampleColor)
                    .traversalDirection(REVERSE)
                    .combinationOperator(OR)
                    .stoppingRule(DestinationStopper.class)
                    .maxBranchLength(1000)
                    .graph(graph)
                    .links(links)
                    .make();

            List<List<List<CortexVertex>>> q = new ArrayList<>();
            for (int i = 0; i < ms.size() - 1; i++) {
                //q.add(Collections.singletonList(ms.get(i)));

                CortexVertex source = ms.get(i).get(ms.get(i).size() - 1);
                CortexVertex sink   = ms.get(i+1).get(0);

                DirectedWeightedPseudograph<CortexVertex, CortexEdge> gFill = ef.dfs(source.getKmerAsString(), sink.getKmerAsString());
                if (gFill == null) {
                    gFill = er.dfs(sink.getKmerAsString(), source.getKmerAsString());
                }

                List<List<CortexVertex>> llc = new ArrayList<>();
                if (gFill != null) {
                    PathFinder pf = new PathFinder(gFill, cs.iterator().next());
                    List<GraphPath<CortexVertex, CortexEdge>> gps = pf.getPaths(source, sink);

                    for (GraphPath<CortexVertex, CortexEdge> gp : gps) {
                        llc.add(gp.getVertexList());
                    }
                }

                if (llc.size() == 0) {
                    q.add(Collections.singletonList(ms.get(i)));
                } else if (llc.size() == 1) {
                    ms.get(i).addAll(llc.get(0));
                    q.add(Collections.singletonList(ms.get(i)));
                } else {
                    q.add(Collections.singletonList(ms.get(i)));
                    q.add(llc);
                }
            }

            q.add(Collections.singletonList(ms.get(ms.size() - 1)));

            mss.put(cs.iterator().next(), q);
        }

        return mss;
    }
    */

    public static DirectedWeightedPseudograph<CortexVertex, CortexEdge> fillGaps(List<CortexVertex> w, CortexGraph graph, List<CortexLinks> links, Set<Integer> colors) {
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> gAll = new DirectedWeightedPseudograph<>(CortexEdge.class);

        for (int c : colors) {
            DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = new DirectedWeightedPseudograph<>(CortexEdge.class);
            for (int i = 0; i < w.size(); i++) {
                CortexVertex v1 = w.get(i);

                if (v1.getCortexRecord().getCoverage(c) > 0) {
                    g.addVertex(v1);

                    //Set<String> anks = convertToStrings(TraversalUtils.getAllPrevKmers(v0.getCortexRecord(), !v0.getKmerAsString().equals(v0.getCanonicalKmer().getKmerAsString())).get(c));
                    if (i > 0 && w.get(i - 1).getCortexRecord().getCoverage(c) > 0) {
                        g.addEdge(w.get(i - 1), v1, new CortexEdge(w.get(i - 1), v1, c, 1.0));
                    }
                }
            }

            Set<String> sources = new HashSet<>();
            Set<String> sinks = new HashSet<>();

            for (CortexVertex v : g.vertexSet()) {
                Set<String> vs = new HashSet<>();
                for (CortexEdge e : g.outgoingEdgesOf(v)) {
                    if (e.getColor() == c) {
                        vs.add(g.getEdgeTarget(e).getKmerAsString());
                    }
                }

                Set<String> remainingNext = convertToStrings(TraversalUtils.getAllNextKmers(v.getCortexRecord(), !v.getKmerAsString().equals(v.getCanonicalKmer().getKmerAsString())).get(c));
                remainingNext.removeAll(vs);

                if (remainingNext.size() > 0) {
                    sources.add(v.getKmerAsString());
                }

                Set<String> vp = new HashSet<>();
                for (CortexEdge e : g.incomingEdgesOf(v)) {
                    if (e.getColor() == c) {
                        vp.add(g.getEdgeSource(e).getKmerAsString());
                    }
                }

                Set<String> remainingPrev = convertToStrings(TraversalUtils.getAllPrevKmers(v.getCortexRecord(), !v.getKmerAsString().equals(v.getCanonicalKmer().getKmerAsString())).get(c));
                remainingPrev.removeAll(vp);

                if (remainingPrev.size() > 0) {
                    sinks.add(v.getKmerAsString());
                }
            }

            TraversalEngine ef = new TraversalEngineFactory()
                    .traversalColors(c)
                    .traversalDirection(FORWARD)
                    .combinationOperator(OR)
                    .stoppingRule(DestinationStopper.class)
                    .maxBranchLength(1000)
                    .graph(graph)
                    .links(links)
                    .make();

            TraversalEngine er = new TraversalEngineFactory()
                    .traversalColors(c)
                    .traversalDirection(REVERSE)
                    .combinationOperator(OR)
                    .stoppingRule(DestinationStopper.class)
                    .maxBranchLength(1000)
                    .graph(graph)
                    .links(links)
                    .make();

            DirectedWeightedPseudograph<CortexVertex, CortexEdge> gFill = ef.dfs(sources, sinks);
            if (gFill == null) {
                gFill = er.dfs(sinks, sources);
            }

            if (gFill != null) {
                Graphs.addGraph(g, gFill);
            }

            Graphs.addGraph(gAll, g);
        }

        return gAll;
    }

    public static DirectedWeightedPseudograph<CortexVertex, CortexEdge> fillGaps(DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, CortexGraph graph, List<CortexLinks> links, Set<Integer> colors) {
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> gFilled = new DirectedWeightedPseudograph<>(CortexEdge.class);
        Graphs.addGraph(gFilled, g);

        Set<String> availableVertices = new HashSet<>();

        for (CortexEdge e : g.edgeSet()) {
            CortexVertex v0 = g.getEdgeSource(e);
            CortexVertex v1 = g.getEdgeTarget(e);

            availableVertices.add(v0.getKmerAsString());
            availableVertices.add(v1.getKmerAsString());

            Map<Integer, Set<CortexByteKmer>> anks = TraversalUtils.getAllNextKmers(v0.getCortexRecord(), !v0.getKmerAsString().equals(v0.getCanonicalKmer().getKmerAsString()));
            for (int c : colors) {
                Set<String> nks = convertToStrings(anks.get(c));

                if (nks.contains(v1.getKmerAsString())) {
                    gFilled.addEdge(v0, v1, new CortexEdge(v0, v1, c, 1.0));
                }
            }
        }

        for (int c : colors) {
            Set<String> sources = new HashSet<>();
            Set<String> sinks = new HashSet<>();

            for (CortexVertex v : gFilled.vertexSet()) {
                Set<String> vs = new HashSet<>();
                for (CortexEdge e : gFilled.outgoingEdgesOf(v)) {
                    if (e.getColor() == c) {
                        vs.add(gFilled.getEdgeTarget(e).getKmerAsString());
                    }
                }

                Set<String> remainingNext = convertToStrings(TraversalUtils.getAllNextKmers(v.getCortexRecord(), !v.getKmerAsString().equals(v.getCanonicalKmer().getKmerAsString())).get(c));
                remainingNext.removeAll(vs);

                if (remainingNext.size() > 0) {
//                    boolean presentInGraph = false;
//                    for (String rn : remainingNext) {
//                        if (availableVertices.contains(rn)) {
//                            presentInGraph = true;
//                            break;
//                        }
//                    }
//
//                    if (presentInGraph) {
                        sources.add(v.getKmerAsString());
//                    }
                }

                Set<String> vp = new HashSet<>();
                for (CortexEdge e : gFilled.incomingEdgesOf(v)) {
                    if (e.getColor() == c) {
                        vp.add(gFilled.getEdgeSource(e).getKmerAsString());
                    }
                }

                Set<String> remainingPrev = convertToStrings(TraversalUtils.getAllPrevKmers(v.getCortexRecord(), !v.getKmerAsString().equals(v.getCanonicalKmer().getKmerAsString())).get(c));
                remainingPrev.removeAll(vp);

                if (remainingPrev.size() > 0) {
//                    boolean presentInGraph = false;
//                    for (String rp : remainingPrev) {
//                        if (availableVertices.contains(rp)) {
//                            presentInGraph = true;
//                            break;
//                        }
//                    }
//
//                    if (presentInGraph) {
                        sinks.add(v.getKmerAsString());
//                    }
                }
            }

            TraversalEngine ef = new TraversalEngineFactory()
                    .traversalColors(c)
                    .traversalDirection(FORWARD)
                    .combinationOperator(OR)
                    .stoppingRule(DestinationStopper.class)
                    .maxBranchLength(1000)
                    .graph(graph)
                    .links(links)
                    .make();

            TraversalEngine er = new TraversalEngineFactory()
                    .traversalColors(c)
                    .traversalDirection(REVERSE)
                    .combinationOperator(OR)
                    .stoppingRule(DestinationStopper.class)
                    .maxBranchLength(1000)
                    .graph(graph)
                    .links(links)
                    .make();

            DirectedWeightedPseudograph<CortexVertex, CortexEdge> gFill = ef.dfs(sources, sinks);
            if (gFill == null) {
                gFill = er.dfs(sinks, sources);
            }

            if (gFill != null) {
                Graphs.addGraph(gFilled, gFill);
            }
        }

        return gFilled;
    }

    private static Set<String> convertToStrings(Set<CortexByteKmer> bs) {
        Set<String> es = new HashSet<>();

        for (CortexByteKmer b : bs) {
            es.add(new String(b.getKmer()));
        }

        return es;
    }

    public static DirectedWeightedPseudograph<CortexVertex, CortexEdge> toGraph(List<CortexVertex> walk, Set<Integer> colors) {
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = new DirectedWeightedPseudograph<>(CortexEdge.class);

        CortexVertex pv = walk.get(0);
        g.addVertex(pv);

        for (int i = 1; i < walk.size(); i++) {
            CortexVertex nv = walk.get(i);

            g.addVertex(nv);

            for (int c : colors) {
                if (pv.getCortexRecord().getCoverage(c) > 0 && nv.getCortexRecord().getCoverage(c) > 0) {
                    g.addEdge(pv, nv, new CortexEdge(pv, nv, c, 1.0));
                }
            }

            pv = nv;
        }

        return g;
    }

    public static DirectedWeightedPseudograph<CortexVertex, CortexEdge> subsetGraph(DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, int color) {
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> gs = new DirectedWeightedPseudograph<>(CortexEdge.class);

        for (CortexEdge e : g.edgeSet()) {
            if (e.getColor() == color) {
                CortexVertex s = g.getEdgeSource(e);
                CortexVertex t = g.getEdgeTarget(e);

                gs.addVertex(s);
                gs.addVertex(t);
                gs.addEdge(s, t, e);
            }
        }

        return gs;
    }

    public static String toContig(List<CortexVertex> walk) {
        StringBuilder sb = new StringBuilder();

        for (CortexVertex cv : walk) {
            String sk = cv.getKmerAsString();

            if (sb.length() == 0) {
                sb.append(sk);
            } else {
                sb.append(sk.substring(sk.length() - 1, sk.length()));
            }
        }

        return sb.toString();
    }

    public static List<CortexVertex> toWalk(DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, CanonicalKmer ck, int color) {
        return TraversalUtils.toWalk(g, ck.getKmerAsString(), color);
    }

    public static List<CortexVertex> toWalk(DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, String sk, int color) {
        List<CortexVertex> w = new ArrayList<>();

        if (g == null) { return w; }

        CortexVertex seed = null;
        for (CortexVertex v : g.vertexSet()) {
            if (v.getKmerAsString().equals(sk) && v.getCortexRecord().getCoverage(color) > 0 && (seed == null || v.getCopyIndex() < seed.getCopyIndex())) {
                seed = v;
            }
        }

        if (seed != null) {
            w.add(seed);

            Set<CortexVertex> seen = new HashSet<>();
            CortexVertex cv = seed;
            while (cv != null && !seen.contains(cv)) {
                List<CortexVertex> nvs = new ArrayList<>();
                for (CortexEdge e : g.outgoingEdgesOf(cv)) {
                    if (e.getColor() == color) {
                        nvs.add(g.getEdgeTarget(e));
                    }
                }

                nvs.remove(cv);

                CortexVertex nv = null;

                if (nvs.size() == 1) {
                    nv = nvs.get(0);
                } else if (nvs.size() > 1) {
                    boolean allKmersTheSame = true;
                    for (int i = 1; i < nvs.size(); i++) {
                        if (!nvs.get(0).getCanonicalKmer().equals(nvs.get(i).getCanonicalKmer())) {
                            allKmersTheSame = false;
                            break;
                        }
                    }

                    if (allKmersTheSame) {
                        nvs.sort((o1, o2) -> o1.getCopyIndex() < o2.getCopyIndex() ? -1 : 1);

                        if (nvs.size() > 0) {
                            nv = nvs.get(0);
                        }
                    }
                }

                if (nv != null) {
                    w.add(nv);
                    seen.add(cv);
                }

                cv = nv;
            }

            seen = new HashSet<>();
            cv = seed;
            while (cv != null && !seen.contains(cv)) {
                List<CortexVertex> pvs = new ArrayList<>();
                for (CortexEdge e : g.incomingEdgesOf(cv)) {
                    if (e.getColor() == color) {
                        pvs.add(g.getEdgeSource(e));
                    }
                }

                pvs.remove(cv);

                CortexVertex pv = null;

                if (pvs.size() == 1) {
                    pv = pvs.get(0);
                } else if (pvs.size() > 1) {
                    boolean allKmersTheSame = true;
                    for (int i = 1; i < pvs.size(); i++) {
                        if (!pvs.get(0).getCanonicalKmer().equals(pvs.get(i).getCanonicalKmer())) {
                            allKmersTheSame = false;
                            break;
                        }
                    }

                    if (allKmersTheSame) {
                        pvs.sort((o1, o2) -> o1.getCopyIndex() > o2.getCopyIndex() ? -1 : 1);

                        if (pvs.size() > 0) {
                            pv = pvs.get(0);
                        }
                    }
                }

                if (pv != null) {
                    w.add(0, pv);
                    seen.add(cv);
                }

                cv = pv;
            }
        }

        return w;
    }

    public static CortexVertex findVertex(DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, CanonicalKmer ck) {
        for (CortexVertex v : g.vertexSet()) {
            if (v.getCanonicalKmer().equals(ck)) {
                return v;
            }
        }

        return null;
    }

    public static CortexVertex findVertex(DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, String sk) {
        for (CortexVertex v : g.vertexSet()) {
            if (v.getKmerAsString().equals(sk)) {
                return v;
            }
        }

        return null;
    }

    public static Map<Integer, Set<CortexByteKmer>> getAllPrevKmers(CortexRecord cr, boolean isFlipped) {
        Map<Integer, Set<CortexByteKmer>> prevKmers = new HashMap<>();

        if (cr != null) {
            byte[] sk = !isFlipped ? cr.getKmerAsBytes() : SequenceUtils.reverseComplement(cr.getKmerAsBytes());
            Map<Integer, Set<Byte>> inEdges = getInEdges(cr, isFlipped);

            for (int c = 0; c < cr.getNumColors(); c++) {
                Set<CortexByteKmer> inKmers = new HashSet<>();

                for (byte inEdge : inEdges.get(c)) {
                    byte[] inKmer = new byte[sk.length];
                    inKmer[0] = inEdge;
                    System.arraycopy(sk, 0, inKmer, 1, sk.length - 1);

                    inKmers.add(new CortexByteKmer(inKmer));
                }

                prevKmers.put(c, inKmers);
            }
        }

        return prevKmers;
    }

    public static Map<Integer, Set<CortexByteKmer>> getAllNextKmers(CortexRecord cr, boolean isFlipped) {
        Map<Integer, Set<CortexByteKmer>> nextKmers = new HashMap<>();

        if (cr != null) {
            byte[] sk = !isFlipped ? cr.getKmerAsBytes() : SequenceUtils.reverseComplement(cr.getKmerAsBytes());
            Map<Integer, Set<Byte>> outEdges = getOutEdges(cr, isFlipped);

            for (int c = 0; c < cr.getNumColors(); c++) {
                Set<CortexByteKmer> outKmers = new HashSet<>();

                for (byte outEdge : outEdges.get(c)) {
                    byte[] outKmer = new byte[sk.length];
                    System.arraycopy(sk, 1, outKmer, 0, sk.length - 1);
                    outKmer[outKmer.length - 1] = outEdge;

                    outKmers.add(new CortexByteKmer(outKmer));
                }

                nextKmers.put(c, outKmers);
            }
        }

        return nextKmers;
    }

    public static Map<Integer, Set<Byte>> getInEdges(CortexRecord cr, boolean kmerIsFlipped) {
        Map<Integer, Set<Byte>> inEdges = new HashMap<>();

        if (!kmerIsFlipped) {
            for (int c = 0; c < cr.getNumColors(); c++) {
                inEdges.put(c, new HashSet<>(cr.getInEdgesAsBytes(c, false)));
            }
        } else {
            for (int c = 0; c < cr.getNumColors(); c++) {
                inEdges.put(c, new HashSet<>(cr.getOutEdgesAsBytes(c, true)));
            }
        }

        return inEdges;
    }

    public static Map<Integer, Set<Byte>> getOutEdges(CortexRecord cr, boolean kmerIsFlipped) {
        Map<Integer, Set<Byte>> outEdges = new HashMap<>();

        if (!kmerIsFlipped) {
            for (int c = 0; c < cr.getNumColors(); c++) {
                outEdges.put(c, new HashSet<>(cr.getOutEdgesAsBytes(c, false)));
            }
        } else {
            for (int c = 0; c < cr.getNumColors(); c++) {
                outEdges.put(c, new HashSet<>(cr.getInEdgesAsBytes(c, true)));
            }
        }

        return outEdges;
    }

    public static int outDegree(DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, CortexVertex v) {
        Set<CortexVertex> vs = new HashSet<>();

        Set<CortexEdge> es = g.outgoingEdgesOf(v);
        for (CortexEdge e : es) {
            vs.add(g.getEdgeTarget(e));
        }

        return vs.size();
    }

    public static int inDegree(DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, CortexVertex v) {
        Set<CortexVertex> vs = new HashSet<>();

        Set<CortexEdge> es = g.incomingEdgesOf(v);
        for (CortexEdge e : es) {
            vs.add(g.getEdgeSource(e));
        }

        return vs.size();
    }
}
