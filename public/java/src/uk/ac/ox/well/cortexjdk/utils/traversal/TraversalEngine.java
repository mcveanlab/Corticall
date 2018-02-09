package uk.ac.ox.well.cortexjdk.utils.traversal;

import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.Nullable;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.Main;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.graph.ConnectivityAnnotations;
import uk.ac.ox.well.cortexjdk.utils.io.graph.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.TraversalStoppingRule;

import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.*;

public class TraversalEngine {
    final private TraversalEngineConfiguration ec;

    public TraversalEngine(TraversalEngineConfiguration ec) { this.ec = ec; }

    private CortexByteKmer curKmer = null;
    private CortexByteKmer prevKmer;
    private CortexByteKmer nextKmer;
    private Set<CortexByteKmer> seen;

    private Set<String> kmerSources;
    private Set<ConnectivityAnnotations> specificLinksFiles;
    private LinkStore linkStore;
    private boolean goForward;

    public TraversalEngineConfiguration getConfiguration() { return ec; }

    private void validateConfiguration(String seed) {
        if (seed.length() != ec.getGraph().getKmerSize()) {
            throw new CortexJDKException("Graph traversal starting kmer is not equal to graph kmer size (" + seed.length() + " vs " + ec.getGraph().getKmerSize() + ")");
        }

        if (ec.getTraversalColor() == -1 || ec.getTraversalColor() >= ec.getGraph().getNumColors()) {
            throw new CortexJDKException("Traversal color '" + ec.getTraversalColor() + "' is invalid.");
        }
    }

    public DirectedWeightedPseudograph<CortexVertex, CortexEdge> dfs(Collection<String> sources) {
        return dfs(sources, null);
    }

    public DirectedWeightedPseudograph<CortexVertex, CortexEdge> dfs(Collection<String> sources, Collection<String> sinks) {
        return null;
    }

    public DirectedWeightedPseudograph<CortexVertex, CortexEdge> dfs(CanonicalKmer source) {
        return dfs(source.getKmerAsString());
    }

    public DirectedWeightedPseudograph<CortexVertex, CortexEdge> dfs(String source) {
        validateConfiguration(source);

        CortexVertex cv = new CortexVertexFactory()
                .bases(source)
                .record(ec.getGraph().findRecord(source))
                .copyIndex(0)
                .make();

        DirectedWeightedPseudograph<CortexVertex, CortexEdge> dfsr = (ec.getTraversalDirection() == BOTH || ec.getTraversalDirection() == REVERSE) ? dfs(cv, false, 0, new HashSet<>()) : null;
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> dfsf = (ec.getTraversalDirection() == BOTH || ec.getTraversalDirection() == FORWARD) ? dfs(cv, true,  0, new HashSet<>()) : null;

        DirectedWeightedPseudograph<CortexVertex, CortexEdge> dfs = null;

        if (ec.getGraphCombinationOperator() == OR) {
            if (dfsr != null || dfsf != null) {
                dfs = new DirectedWeightedPseudograph<>(CortexEdge.class);

                if (dfsr != null) { Graphs.addGraph(dfs, dfsr); }
                if (dfsf != null) { Graphs.addGraph(dfs, dfsf); }
            }
        } else {
            if (dfsr != null && dfsf != null) {
                dfs = new DirectedWeightedPseudograph<>(CortexEdge.class);

                Graphs.addGraph(dfs, dfsr);
                Graphs.addGraph(dfs, dfsf);
            }
        }

        if (dfs != null) {
            return addSecondaryColors(dfs);
        }

        return null;
    }

    public List<CortexVertex> walk(CanonicalKmer seed) { return toWalk(dfs(seed.getKmerAsString()), seed.getKmerAsString()); }

    public List<CortexVertex> walk(String seed) {
        return toWalk(dfs(seed), seed);
    }

    public List<CortexVertex> assemble(String seed) {
        List<CortexVertex> contig = new ArrayList<>();
        CortexVertex sv = new CortexVertex(new CortexByteKmer(seed.getBytes()), ec.getGraph().findRecord(seed));
        contig.add(sv);

        contig.addAll(assemble(seed, true));
        contig.addAll(0, assemble(seed, false));

        return contig;
    }

    public List<CortexVertex> assemble(String seed, boolean goForward) {
        List<CortexVertex> contig = new ArrayList<>();

        seek(seed);
        if (goForward) {
            while (hasNext() && contig.size() < ec.getMaxBranchLength()) {
                CortexVertex cv = next();
                contig.add(cv);
            }
        } else {
            while (hasPrevious() && contig.size() < ec.getMaxBranchLength()) {
                CortexVertex cv = previous();
                contig.add(0, cv);
            }
        }

        return contig;
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

    public static List<CortexVertex> toWalk(DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, CanonicalKmer ck) {
        return toWalk(g, ck.getKmerAsString());
    }

    public static List<CortexVertex> toWalk(DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, String sk) {
        List<CortexVertex> w = new ArrayList<>();

        if (g == null) { return w; }

        CortexVertex seed = null;
        for (CortexVertex v : g.vertexSet()) {
            if (v.getKmerAsString().equals(sk) && v.getCopyIndex() == 0) {
                seed = v;
                break;
            }
        }

        if (seed != null) {
            w.add(seed);

            Set<CortexVertex> seen = new HashSet<>();
            CortexVertex cv = seed;
            while (cv != null && !seen.contains(cv)) {
                List<CortexVertex> nvs = Graphs.successorListOf(g, cv);
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
                List<CortexVertex> pvs = Graphs.predecessorListOf(g, cv);
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

    private static TraversalStoppingRule<CortexVertex, CortexEdge> instantiateStopper(Class<? extends TraversalStoppingRule<CortexVertex, CortexEdge>> stopperClass) {
        try {
            return stopperClass.newInstance();
        } catch (InstantiationException e) {
            throw new CortexJDKException("Could not instantiate stoppingRule: ", e);
        } catch (IllegalAccessException e) {
            throw new CortexJDKException("Illegal access while trying to instantiate stoppingRule: ", e);
        }
    }

    @Nullable
    private DirectedWeightedPseudograph<CortexVertex, CortexEdge> dfs(CortexVertex cv, boolean goForward, int currentJunctionDepth, Set<CortexVertex> visitedOld) {
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = new DirectedWeightedPseudograph<>(CortexEdge.class);

        // Account for vertices visited in progenitor branches (but not other progeny branches)
        Set<CortexVertex> visited = new HashSet<>(visitedOld);

        // If links are available, reset the state of the LinkStore
        if (!ec.getLinks().isEmpty()) {
            seek(cv.getKmerAsString());
        }

        // Instantiate a new stopping rule per branch
        TraversalStoppingRule<CortexVertex, CortexEdge> stoppingRule = instantiateStopper(ec.getStoppingRule());

        Set<CortexVertex> avs;

        do {
            Set<CortexVertex> pvs = getPrevVertices(cv.getKmerAsByteKmer());
            Set<CortexVertex> nvs = getNextVertices(cv.getKmerAsByteKmer());
            avs = goForward ? nvs : pvs;

            if (!ec.getLinks().isEmpty()) {
                // If we have links, then we are permitted to traverse some vertices multiple times.  Include a copy
                // count for those vertices so that we can distinguish each copy in the resulting subgraph.

                CortexVertex qv = null;
                if (goForward && hasNext()) {
                    qv = next();
                } else if (!goForward && hasPrevious()) {
                    qv = previous();
                }

                if (qv != null) {
                    CortexVertex lv = null;
                    do {
                        int copyIndex;
                        if (goForward) { copyIndex = lv == null ? 1 : lv.getCopyIndex() + 1; }
                        else { copyIndex = lv == null ? -1 : lv.getCopyIndex() - 1; }

                        lv = new CortexVertexFactory()
                                .bases(qv.getKmerAsString())
                                .record(qv.getCortexRecord())
                                .copyIndex(copyIndex)
                                .make();
                    } while (visited.contains(lv));

                    avs = new HashSet<>();
                    avs.add(lv);
                }
            }

            // Connect all neighboring vertices to the graph (useful for visualization)
            if (ec.connectAllNeighbors()) {
                connectVertex(g, cv, pvs, nvs);
            }

            // Avoid traversing infinite loops by removing from traversal consideration
            // those vertices that have already been incorporated into the graph.
            Set<CortexVertex> seen = new HashSet<>();
            for (CortexVertex av : avs) {
                if (visited.contains(av)) {
                    seen.add(av);
                }
            }
            avs.removeAll(seen);

            boolean previouslyVisited = visited.contains(cv);
            visited.add(cv);

            // Decide if we should keep exploring the graph or not
            TraversalState<CortexVertex> ts = new TraversalState<>(cv, goForward, ec.getTraversalColor(), ec.getJoiningColors(), currentJunctionDepth, g.vertexSet().size(), avs.size(), false, g.vertexSet().size() > ec.getMaxBranchLength(), ec.getSink(), ec.getRois());

            if (!previouslyVisited && stoppingRule.keepGoing(ts)) {
                if (avs.size() == 1) {
                    if (goForward) { connectVertex(g, cv, null, avs);  }
                    else           { connectVertex(g, cv, avs,  null); }

                    cv = avs.iterator().next();
                } else {
                    boolean childrenWereSuccessful = false;

                    for (CortexVertex av : avs) {
                        DirectedWeightedPseudograph<CortexVertex, CortexEdge> branch = dfs(av, goForward, currentJunctionDepth + 1, visited);

                        if (branch != null) {
                            if (goForward) { connectVertex(branch, cv, null, Collections.singleton(av));  }
                            else           { connectVertex(branch, cv, Collections.singleton(av),  null); }

                            Graphs.addGraph(g, branch);
                            childrenWereSuccessful = true;
                        } else {
                            // could mark a rejected traversal here rather than just throwing it away
                        }
                    }

                    TraversalState<CortexVertex> tsChild = new TraversalState<>(cv, goForward, ec.getTraversalColor(), ec.getJoiningColors(), currentJunctionDepth, g.vertexSet().size(), avs.size(), true, g.vertexSet().size() > ec.getMaxBranchLength(), ec.getSink(), ec.getRois());

                    if (childrenWereSuccessful || stoppingRule.hasTraversalSucceeded(tsChild)) {
                        return g;
                    } else {
                        // could mark a rejected traversal here rather than just throwing it away
                    }
                }
            } else if (stoppingRule.traversalSucceeded()) {
                return g;
            } else {
                return null;
            }
        } while (avs.size() == 1);

        return null;
    }

    private void connectVertex(DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, CortexVertex cv, Set<CortexVertex> pvs, Set<CortexVertex> nvs) {
        g.addVertex(cv);

        if (pvs != null) {
            for (CortexVertex pv : pvs) {
                g.addVertex(pv);

                if (!g.containsEdge(pv, cv)) {
                    g.addEdge(pv, cv, new CortexEdge(ec.getTraversalColor(), 1.0));
                }
            }
        }

        if (nvs != null) {
            for (CortexVertex nv : nvs) {
                g.addVertex(nv);

                if (!g.containsEdge(cv, nv)) {
                    g.addEdge(cv, nv, new CortexEdge(ec.getTraversalColor(), 1.0));
                }
            }
        }
    }

    public Set<CortexVertex> getPrevVertices(CortexByteKmer sk) {
        Set<CortexVertex> prevVertices = new HashSet<>();

        Map<Integer, Set<CortexByteKmer>> prevKmers = getAllPrevKmers(sk);

        if (prevKmers.get(ec.getTraversalColor()).size() > 0) {
            for (CortexByteKmer prevKmer : prevKmers.get(ec.getTraversalColor())) {
                prevVertices.add(new CortexVertex(prevKmer, ec.getGraph().findRecord(prevKmer)));
            }
        } else {
            Map<CortexByteKmer, Set<Integer>> inKmerMap = new HashMap<>();

            for (int c : ec.getRecruitmentColors()) {
                Set<CortexByteKmer> inKmers = prevKmers.get(c);

                for (CortexByteKmer inKmer : inKmers) {
                    if (!inKmerMap.containsKey(inKmer)) {
                        inKmerMap.put(inKmer, new HashSet<>());
                    }

                    inKmerMap.get(inKmer).add(c);
                }
            }

            for (CortexByteKmer prevKmer : inKmerMap.keySet()) {
                prevVertices.add(new CortexVertex(prevKmer, ec.getGraph().findRecord(prevKmer)));
            }
        }

        return prevVertices;
    }

    public Set<CortexVertex> getNextVertices(CortexByteKmer sk) {
        Set<CortexVertex> nextVertices = new HashSet<>();

        Map<Integer, Set<CortexByteKmer>> nextKmers = getAllNextKmers(sk);
        if (nextKmers.get(ec.getTraversalColor()).size() > 0) {
            for (CortexByteKmer nextKmer : nextKmers.get(ec.getTraversalColor())) {
                nextVertices.add(new CortexVertex(nextKmer, ec.getGraph().findRecord(nextKmer)));
            }
        } else {
            Map<CortexByteKmer, Set<Integer>> outKmerMap = new HashMap<>();

            for (int c : ec.getRecruitmentColors()) {
                Set<CortexByteKmer> outKmers = nextKmers.get(c);

                for (CortexByteKmer outKmer : outKmers) {
                    if (!outKmerMap.containsKey(outKmer)) {
                        outKmerMap.put(outKmer, new HashSet<>());
                    }

                    outKmerMap.get(outKmer).add(c);
                }
            }

            for (CortexByteKmer nextKmer : outKmerMap.keySet()) {
                nextVertices.add(new CortexVertex(nextKmer, ec.getGraph().findRecord(nextKmer)));
            }
        }

        return nextVertices;
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

    public Map<Integer, Set<CortexByteKmer>> getAllPrevKmers(CortexByteKmer sk) {
        CanonicalKmer ck = new CanonicalKmer(sk.getKmer());
        CortexRecord cr = ec.getGraph().findRecord(ck);

        return getAllPrevKmers(cr, ck.isFlipped());
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

    public Map<Integer, Set<CortexByteKmer>> getAllNextKmers(CortexByteKmer sk) {
        CanonicalKmer ck = new CanonicalKmer(sk.getKmer());
        CortexRecord cr = ec.getGraph().findRecord(ck);

        return getAllNextKmers(cr, ck.isFlipped());
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

    public CortexVertex next() {
        if (nextKmer == null) { throw new NoSuchElementException("No single advance kmer from cursor '" + curKmer + "'"); }
        if (specificLinksFiles == null || !goForward) {
            goForward = true;
            seek(new String(curKmer.getKmer()));

            initializeLinkStore(goForward);
        }

        updateLinkStore(goForward);

        CortexRecord cr = ec.getGraph().findRecord(nextKmer);
        CortexVertex cv = new CortexVertex(nextKmer, cr, kmerSources);

        prevKmer = curKmer;
        curKmer = nextKmer;

        Set<CortexVertex> nextKmers = getNextVertices(curKmer);
        nextKmer = null;
        kmerSources = null;

        if (nextKmers.size() == 1 && (!seen.contains(nextKmers.iterator().next().getKmerAsByteKmer()) || linkStore.isActive())) {
            nextKmer = nextKmers.iterator().next().getKmerAsByteKmer();

            seen.add(nextKmer);
        } else if (nextKmers.size() > 1) {
            Pair<CortexByteKmer, Set<String>> akp = getAdjacentKmer(curKmer, nextKmers, goForward);
            nextKmer = akp != null ? akp.getFirst() : null;
            kmerSources = akp != null ? akp.getSecond() : null;

            linkStore.incrementAges();
        }

        if (linkStore.numNewPaths() > 0) {
            linkStore.incrementAges();
        }

        return cv;
    }

    public CortexVertex previous() {
        if (prevKmer == null) { throw new NoSuchElementException("No single prev kmer from cursor '" + curKmer + "'"); }
        if (specificLinksFiles == null || goForward) {
            goForward = false;
            seek(new String(curKmer.getKmer()));

            initializeLinkStore(goForward);
        }

        updateLinkStore(goForward);

        CortexRecord cr = ec.getGraph().findRecord(prevKmer);
        CortexVertex cv = new CortexVertex(prevKmer, cr, kmerSources);

        nextKmer = curKmer;
        curKmer = prevKmer;

        Set<CortexVertex> prevKmers = getPrevVertices(curKmer);
        prevKmer = null;
        kmerSources = null;

        if (prevKmers.size() == 1 && (!seen.contains(prevKmers.iterator().next().getKmerAsByteKmer()) || linkStore.isActive())) {
            prevKmer = prevKmers.iterator().next().getKmerAsByteKmer();

            seen.add(prevKmer);
        } else if (prevKmers.size() > 1) {
            Pair<CortexByteKmer, Set<String>> akp = getAdjacentKmer(curKmer, prevKmers, goForward);
            prevKmer = akp != null ? akp.getFirst() : null;
            kmerSources = akp != null ? akp.getSecond() : null;

            linkStore.incrementAges();
        }

        if (linkStore.numNewPaths() > 0) {
            linkStore.incrementAges();
        }

        return cv;
    }

    private Pair<CortexByteKmer, Set<String>> getAdjacentKmer(CortexByteKmer kmer, Set<CortexVertex> adjKmers, boolean goForward) {
        Pair<String, Set<String>> choicePair = linkStore.getNextJunctionChoice();
        String choice = choicePair.getFirst();
        Set<String> sources = choicePair.getSecond();

        if (choice != null) {
            Set<CortexByteKmer> adjKmersStr = new HashSet<>();
            for (CortexVertex cv : adjKmers) {
                adjKmersStr.add(cv.getKmerAsByteKmer());
            }

            byte[] bAdjKmer = new byte[kmer.length()];
            if (goForward) {
                System.arraycopy(kmer.getKmer(), 1, bAdjKmer, 0, kmer.length() - 1);
                bAdjKmer[bAdjKmer.length - 1] = choice.getBytes()[0];
            } else {
                bAdjKmer[0] = choice.getBytes()[0];
                System.arraycopy(kmer.getKmer(), 0, bAdjKmer, 1, kmer.length() - 1);
            }

            CortexByteKmer adjKmer = new CortexByteKmer(bAdjKmer);

            if (adjKmersStr.contains(adjKmer)) {
                return new Pair<>(adjKmer, sources);
            }
        }

        return null;
    }

    private void initializeLinkStore(boolean goForward) {
        specificLinksFiles = new HashSet<>();

        if (!ec.getLinks().isEmpty()) {
            for (ConnectivityAnnotations lm : ec.getLinks()) {
                if (lm.getHeader().getSampleNameForColor(0).equals(ec.getGraph().getSampleName(ec.getTraversalColor()))) {
                    specificLinksFiles.add(lm);

                    CanonicalKmer ck = new CanonicalKmer(curKmer.getKmer());
                    if (lm.containsKey(ck)) {
                        linkStore.add(curKmer, lm.get(ck), goForward, lm.getSource());
                    }
                }
            }
        }
    }

    private void updateLinkStore(boolean goForward) {
        specificLinksFiles = new HashSet<>();
        if (!ec.getLinks().isEmpty()) {
            for (ConnectivityAnnotations lm : ec.getLinks()) {
                if (lm.getHeader().getSampleNameForColor(0).equals(ec.getGraph().getSampleName(ec.getTraversalColor()))) {
                    specificLinksFiles.add(lm);

                    if (goForward) {
                        CanonicalKmer nk = nextKmer == null ? null : new CanonicalKmer(nextKmer.getKmer());
                        if (nextKmer != null && lm.containsKey(nk)) {
                            linkStore.add(nextKmer, lm.get(nk), goForward, lm.getSource());
                        }
                    } else {
                        CanonicalKmer pk = prevKmer == null ? null : new CanonicalKmer(prevKmer.getKmer());
                        if (prevKmer != null && lm.containsKey(pk)) {
                            linkStore.add(prevKmer, lm.get(pk), goForward, lm.getSource());
                        }
                    }
                }
            }
        }
    }

    public void seek(String sk) {
        if (sk != null) {
            curKmer = new CortexByteKmer(sk.getBytes());

            Set<CortexVertex> prevKmers = getPrevVertices(curKmer);
            prevKmer = (prevKmers.size() == 1) ? prevKmers.iterator().next().getKmerAsByteKmer() : null;

            Set<CortexVertex> nextKmers = getNextVertices(curKmer);
            nextKmer = (nextKmers.size() == 1) ? nextKmers.iterator().next().getKmerAsByteKmer() : null;

            linkStore = new LinkStore();
            seen = new HashSet<>();
            specificLinksFiles = null;
        }
    }

    public boolean hasNext() { return nextKmer != null; }

    public boolean hasPrevious() { return prevKmer != null; }

    private DirectedWeightedPseudograph<CortexVertex, CortexEdge> addSecondaryColors(DirectedWeightedPseudograph<CortexVertex, CortexEdge> g) {
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> m = new DirectedWeightedPseudograph<>(CortexEdge.class);
        Graphs.addGraph(m, g);

        if (!ec.getSecondaryColors().isEmpty()) {
            Map<CortexByteKmer, Map<Integer, Set<CortexByteKmer>>> pkscache = new HashMap<>();
            Map<CortexByteKmer, Map<Integer, Set<CortexByteKmer>>> nkscache = new HashMap<>();

            for (int c : ec.getSecondaryColors()) {
                if (c != ec.getTraversalColor()) {
                    DirectedWeightedPseudograph<CortexVertex, CortexEdge> g2 = new DirectedWeightedPseudograph<>(CortexEdge.class);

                    for (CortexVertex v : g.vertexSet()) {
                        Map<Integer, Set<CortexByteKmer>> pks = pkscache.containsKey(v.getKmerAsByteKmer()) ? pkscache.get(v.getKmerAsByteKmer()) : getAllPrevKmers(v.getKmerAsByteKmer());
                        Map<Integer, Set<CortexByteKmer>> nks = nkscache.containsKey(v.getKmerAsByteKmer()) ? nkscache.get(v.getKmerAsByteKmer()) : getAllNextKmers(v.getKmerAsByteKmer());

                        pkscache.put(v.getKmerAsByteKmer(), pks);
                        nkscache.put(v.getKmerAsByteKmer(), nks);

                        g2.addVertex(v);

                        for (CortexByteKmer pk : pks.get(c)) {
                            CortexVertex pv = new CortexVertex(pk, ec.getGraph().findRecord(pk));

                            g2.addVertex(pv);
                            if (!g2.containsEdge(pv, v)) {
                                g2.addEdge(pv, v, new CortexEdge(c, 1.0));
                            }
                        }

                        for (CortexByteKmer nk : nks.get(c)) {
                            CortexVertex nv = new CortexVertex(nk, ec.getGraph().findRecord(nk));

                            g2.addVertex(nv);
                            if (!g2.containsEdge(v, nv)) {
                                g2.addEdge(v, nv, new CortexEdge(c, 1.0));
                            }
                        }
                    }

                    Graphs.addGraph(m, g2);
                }
            }
        }

        return m;
    }
}
