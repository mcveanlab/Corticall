package uk.ac.ox.well.indiana.utils.traversal;

import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.Nullable;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.ByteKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;
import uk.ac.ox.well.indiana.utils.stoppingconditions.TraversalStopper;

import java.util.*;

import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.REVERSE;

public class TraversalEngine implements ListIterator<CortexVertex> {
    private TraversalEngineConfiguration ec;

    public TraversalEngine(TraversalEngineConfiguration ec) { this.ec = ec; }

    public TraversalEngineConfiguration getConfiguration() { return ec; }

    private ByteKmer curKmer = null;
    private ByteKmer prevKmer;
    private ByteKmer nextKmer;
    private Set<ByteKmer> seen;

    private Set<String> kmerSources;

    private LinkStore linkStore;

    public DirectedWeightedPseudograph<CortexVertex, CortexEdge> dfs(String sk) {
        if (sk.length() != ec.getGraph().getKmerSize()) {
            throw new IndianaException("Graph traversal starting kmer is not equal to graph kmer size (" + sk.length() + " vs " + ec.getGraph().getKmerSize() + ")");
        }

        DirectedGraph<CortexVertex, CortexEdge> dfsr = (ec.getTraversalDirection() == BOTH || ec.getTraversalDirection() == REVERSE) ? dfs(sk, false, 0, new HashSet<>()) : null;
        DirectedGraph<CortexVertex, CortexEdge> dfsf = (ec.getTraversalDirection() == BOTH || ec.getTraversalDirection() == FORWARD) ? dfs(sk, true,  0, new HashSet<>()) : null;

        DirectedGraph<CortexVertex, CortexEdge> dfs = null;

        if (ec.getGraphCombinationOperator() == OR) {
            if (dfsr != null || dfsf != null) {
                dfs = new DefaultDirectedGraph<>(CortexEdge.class);

                if (dfsr != null) { Graphs.addGraph(dfs, dfsr); }
                if (dfsf != null) { Graphs.addGraph(dfs, dfsf); }
            }
        } else {
            if (dfsr != null && dfsf != null) {
                dfs = new DefaultDirectedGraph<>(CortexEdge.class);

                Graphs.addGraph(dfs, dfsr);
                Graphs.addGraph(dfs, dfsf);
            }
        }

        /*
        if (dfs != null) {
            return addDisplayColors(dfs);
        }
        */

        return null;
    }

    public List<CortexVertex> walk(String seed) {
        if (ec.getTraversalColor() < 0) { throw new IndianaException("Traversal color not set."); }

        List<CortexVertex> contig = new ArrayList<>();

        CortexVertex sv = new CortexVertex(new ByteKmer(seed.getBytes()), ec.getGraph().findRecord(seed));
        contig.add(sv);

        setCursor(seed, true);
        while (hasNext()) {
            CortexVertex cv = next();
            contig.add(cv);
        }

        setCursor(seed, false);
        while (hasPrevious()) {
            CortexVertex cv = previous();
            contig.add(0, cv);
        }

        return contig;
    }

    public static String toContig(List<CortexVertex> walk) {
        StringBuilder sb = new StringBuilder();

        for (CortexVertex cv : walk) {
            String sk = cv.getSk();

            if (sb.length() == 0) {
                sb.append(sk);
            } else {
                sb.append(sk.substring(sk.length() - 1, sk.length()));
            }
        }

        return sb.toString();
    }

    public static DirectedWeightedPseudograph<CortexVertex, CortexEdge> toPseudograph(CortexGraph graph, List<CortexVertex> walk, int altColor, int ... refColors) {
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> dwp = new DirectedWeightedPseudograph<>(CortexEdge.class);

        dwp.addVertex(walk.get(0));

        for (int i = 1; i < walk.size(); i++) {
            CortexVertex cv0 = walk.get(i - 1);
            CortexVertex cv1 = walk.get(i);

            dwp.addVertex(cv1);
            dwp.addEdge(cv0, cv1, new CortexEdge(altColor, 1.0));
        }

        for (CortexVertex cv : walk) {
            Map<Integer, Set<ByteKmer>> nextKmers = TraversalEngine.getAllNextKmers(cv.getCr(), !cv.getCk().getKmerAsString().equals(cv.getSk()));
            Map<Integer, Set<ByteKmer>> prevKmers = TraversalEngine.getAllPrevKmers(cv.getCr(), !cv.getCk().getKmerAsString().equals(cv.getSk()));

            for (int refColor : refColors) {
                for (ByteKmer nextKmer : nextKmers.get(refColor)) {
                    CortexVertex nv = new CortexVertex(nextKmer, graph.findRecord(nextKmer));

                    dwp.addVertex(nv);
                    dwp.addEdge(cv, nv, new CortexEdge(refColor, 1.0));
                }

                for (ByteKmer prevKmer : prevKmers.get(refColor)) {
                    CortexVertex pv = new CortexVertex(prevKmer, graph.findRecord(prevKmer));

                    dwp.addVertex(pv);
                    dwp.addEdge(pv, cv, new CortexEdge(refColor, 1.0));
                }
            }
        }

        return dwp;
    }

    private static TraversalStopper<CortexVertex, CortexEdge> instantiateStopper(Class<? extends TraversalStopper<CortexVertex, CortexEdge>> stopperClass) {
        try {
            return stopperClass.newInstance();
        } catch (InstantiationException e) {
            throw new IndianaException("Could not instantiate stopper: ", e);
        } catch (IllegalAccessException e) {
            throw new IndianaException("Illegal access while trying to instantiate stopper: ", e);
        }
    }

    @Nullable
    private DirectedGraph<CortexVertex, CortexEdge> dfs(String sk, boolean goForward, int currentTraversalDepth, Set<CortexVertex> visited) {
        DirectedGraph<CortexVertex, CortexEdge> g = new DefaultDirectedGraph<>(CortexEdge.class);

        CortexRecord cr = ec.getGraph().findRecord(new CortexKmer(sk));

        if (cr == null) { throw new IndianaException("Record '" + sk + "' does not exist in graph."); }

        CortexVertex cv = new CortexVertex(new ByteKmer(sk), cr);

        Set<CortexVertex> avs;

        TraversalStopper<CortexVertex, CortexEdge> stoppingRule = instantiateStopper(ec.getStoppingRule());

        do {
            Set<CortexVertex> pvs = getPrevVertices(cv.getBk());
            Set<CortexVertex> nvs = getNextVertices(cv.getBk());
            avs = goForward ? nvs : pvs;

            // Connect the new vertices to the graph
            if (ec.connectAllNeighbors()) {
                connectVertex(g, cv, pvs,  nvs);
            } else if (goForward) {
                connectVertex(g, cv, null, nvs);
            } else {
                connectVertex(g, cv, pvs,  null);
            }

            visited.add(cv);

            // Avoid traversing infinite loops by removing from traversal consideration
            // those vertices that have already been incorporated into the graph.
            Set<CortexVertex> seen = new HashSet<>();
            for (CortexVertex av : avs) {
                if (visited.contains(av)) {
                    seen.add(av);
                }
            }
            avs.removeAll(seen);

            // Decide if we should keep exploring the graph or not
            if (stoppingRule.keepGoing(cv, goForward, ec.getTraversalColor(), ec.getJoiningColors(), currentTraversalDepth, g.vertexSet().size(), avs.size(), false, ec.getPreviousTraversal(), ec.getRois())) {
                if (avs.size() == 1) {
                    cv = avs.iterator().next();
                } else if (avs.size() != 1) {
                    boolean childrenWereSuccessful = false;

                    for (CortexVertex av : avs) {
                        DirectedGraph<CortexVertex, CortexEdge> branch = dfs(av.getSk(), goForward, currentTraversalDepth + 1, visited);

                        if (branch != null) {
                            Graphs.addGraph(g, branch);
                            childrenWereSuccessful = true;
                        } else {
                            // could mark a rejected traversal here rather than just throwing it away
                        }
                    }

                    if (childrenWereSuccessful || stoppingRule.hasTraversalSucceeded(cv, goForward, ec.getTraversalColor(), ec.getJoiningColors(), currentTraversalDepth, g.vertexSet().size(), avs.size(), true, ec.getPreviousTraversal(), ec.getRois())) {
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

    private void connectVertex(DirectedGraph<CortexVertex, CortexEdge> g, CortexVertex cv, Set<CortexVertex> pvs, Set<CortexVertex> nvs) {
        g.addVertex(cv);

        if (pvs != null) {
            for (CortexVertex pv : pvs) {
                g.addVertex(pv);
                g.addEdge(pv, cv, new CortexEdge(ec.getTraversalColor(), 1.0));
            }
        }

        if (nvs != null) {
            for (CortexVertex nv : nvs) {
                g.addVertex(nv);
                g.addEdge(cv, nv, new CortexEdge(ec.getTraversalColor(), 1.0));
            }
        }
    }

    public Set<CortexVertex> getPrevVertices(ByteKmer sk) {
        Set<CortexVertex> prevVertices = new HashSet<>();

        Map<Integer, Set<ByteKmer>> prevKmers = getAllPrevKmers(sk);

        if (prevKmers.get(ec.getTraversalColor()).size() > 0) {
            for (ByteKmer prevKmer : prevKmers.get(ec.getTraversalColor())) {
                prevVertices.add(new CortexVertex(prevKmer, ec.getGraph().findRecord(prevKmer)));
            }
        } else {
            Map<ByteKmer, Set<Integer>> inKmerMap = new HashMap<>();

            for (int c : ec.getRecruitmentColors()) {
                Set<ByteKmer> inKmers = prevKmers.get(c);

                for (ByteKmer inKmer : inKmers) {
                    if (!inKmerMap.containsKey(inKmer)) {
                        inKmerMap.put(inKmer, new HashSet<>());
                    }

                    inKmerMap.get(inKmer).add(c);
                }
            }

            for (ByteKmer prevKmer : inKmerMap.keySet()) {
                prevVertices.add(new CortexVertex(prevKmer, ec.getGraph().findRecord(prevKmer)));
            }
        }

        return prevVertices;
    }

    //public int getInDegree(String sk) { return getPrevVertices(sk).size(); }

    //public int getOutDegree(String sk) { return getNextVertices(sk).size(); }

    public Set<CortexVertex> getNextVertices(ByteKmer sk) {
        Set<CortexVertex> nextVertices = new HashSet<>();

        Map<Integer, Set<ByteKmer>> nextKmers = getAllNextKmers(sk);
        if (nextKmers.get(ec.getTraversalColor()).size() > 0) {
            for (ByteKmer nextKmer : nextKmers.get(ec.getTraversalColor())) {
                nextVertices.add(new CortexVertex(nextKmer, ec.getGraph().findRecord(nextKmer)));
            }
        } else {
            Map<ByteKmer, Set<Integer>> outKmerMap = new HashMap<>();

            for (int c : ec.getRecruitmentColors()) {
                Set<ByteKmer> outKmers = nextKmers.get(c);

                for (ByteKmer outKmer : outKmers) {
                    if (!outKmerMap.containsKey(outKmer)) {
                        outKmerMap.put(outKmer, new HashSet<>());
                    }

                    outKmerMap.get(outKmer).add(c);
                }
            }

            for (ByteKmer nextKmer : outKmerMap.keySet()) {
                nextVertices.add(new CortexVertex(nextKmer, ec.getGraph().findRecord(nextKmer)));
            }
        }

        return nextVertices;
    }

    public static Map<Integer, Set<ByteKmer>> getAllPrevKmers(CortexRecord cr, boolean isFlipped) {
        Map<Integer, Set<ByteKmer>> prevKmers = new HashMap<>();

        if (cr != null) {
            byte[] sk = !isFlipped ? cr.getKmerAsBytes() : SequenceUtils.reverseComplement(cr.getKmerAsBytes());
            Map<Integer, Set<Byte>> inEdges = getInEdges(cr, isFlipped);

            for (int c = 0; c < cr.getNumColors(); c++) {
                Set<ByteKmer> inKmers = new HashSet<>();

                for (byte inEdge : inEdges.get(c)) {
                    byte[] inKmer = new byte[sk.length];
                    inKmer[0] = inEdge;
                    System.arraycopy(sk, 0, inKmer, 1, sk.length - 1);

                    inKmers.add(new ByteKmer(inKmer));
                }

                prevKmers.put(c, inKmers);
            }
        }

        return prevKmers;
    }

    public Map<Integer, Set<ByteKmer>> getAllPrevKmers(ByteKmer sk) {
        CortexKmer ck = new CortexKmer(sk.getKmer());
        CortexRecord cr = ec.getGraph().findRecord(ck);

        return getAllPrevKmers(cr, ck.isFlipped());
    }

    public static Map<Integer, Set<ByteKmer>> getAllNextKmers(CortexRecord cr, boolean isFlipped) {
        Map<Integer, Set<ByteKmer>> nextKmers = new HashMap<>();

        if (cr != null) {
            byte[] sk = !isFlipped ? cr.getKmerAsBytes() : SequenceUtils.reverseComplement(cr.getKmerAsBytes());
            Map<Integer, Set<Byte>> outEdges = getOutEdges(cr, isFlipped);

            for (int c = 0; c < cr.getNumColors(); c++) {
                Set<ByteKmer> outKmers = new HashSet<>();

                for (byte outEdge : outEdges.get(c)) {
                    byte[] outKmer = new byte[sk.length];
                    System.arraycopy(sk, 1, outKmer, 0, sk.length - 1);
                    outKmer[outKmer.length - 1] = outEdge;

                    outKmers.add(new ByteKmer(outKmer));
                }

                nextKmers.put(c, outKmers);
            }
        }

        return nextKmers;
    }

    public Map<Integer, Set<ByteKmer>> getAllNextKmers(ByteKmer sk) {
        CortexKmer ck = new CortexKmer(sk.getKmer());
        CortexRecord cr = ec.getGraph().findRecord(ck);

        return getAllNextKmers(cr, ck.isFlipped());
    }

    private static Map<Integer, Set<Byte>> getInEdges(CortexRecord cr, boolean kmerIsFlipped) {
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

    private static Map<Integer, Set<Byte>> getOutEdges(CortexRecord cr, boolean kmerIsFlipped) {
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

    @Override
    public CortexVertex next() {
        if (nextKmer == null) { throw new NoSuchElementException("No single next kmer from cursor '" + curKmer + "'"); }

        linkStore.incrementAge();

        CortexRecord cr = ec.getGraph().findRecord(nextKmer);
        CortexVertex cv = new CortexVertex(nextKmer, cr, kmerSources);

        prevKmer = curKmer;
        curKmer = nextKmer;

        Set<CortexVertex> nextKmers = getNextVertices(curKmer);
        nextKmer = null;
        kmerSources = null;

        if (nextKmers.size() == 1 && (!seen.contains(nextKmers.iterator().next().getBk()) || linkStore.isActive())) {
            nextKmer = nextKmers.iterator().next().getBk();

            seen.add(nextKmer);
        } else if (nextKmers.size() > 1) {
            Pair<ByteKmer, Set<String>> akp = getAdjacentKmer(curKmer, nextKmers, true);
            nextKmer = akp != null ? akp.getFirst() : null;
            kmerSources = akp != null ? akp.getSecond() : null;
        }

        if (!ec.getLinks().isEmpty()) {
            for (CortexLinks lm : ec.getLinks().keySet()) {
                CortexKmer ck = new CortexKmer(curKmer.getKmer());
                if (lm.containsKey(ck)) {
                    linkStore.add(curKmer, lm.get(ck), true, ec.getLinks().get(lm));
                }
            }
        }

        return cv;
    }

    @Override
    public CortexVertex previous() {
        if (prevKmer == null) { throw new NoSuchElementException("No single prev kmer from cursor '" + curKmer + "'"); }

        linkStore.incrementAge();

        CortexRecord cr = ec.getGraph().findRecord(new CortexKmer(prevKmer.getKmer()));
        CortexVertex cv = new CortexVertex(prevKmer, cr, kmerSources);

        nextKmer = curKmer;
        curKmer = prevKmer;

        Set<CortexVertex> prevKmers = getPrevVertices(curKmer);
        prevKmer = null;
        kmerSources = null;

        if (prevKmers.size() == 1 && (!seen.contains(prevKmers.iterator().next().getBk()) || linkStore.isActive())) {
            prevKmer = prevKmers.iterator().next().getBk();

            seen.add(prevKmer);
        } else if (prevKmers.size() > 1) {
            Pair<ByteKmer, Set<String>> akp = getAdjacentKmer(curKmer, prevKmers, false);
            prevKmer = akp != null ? akp.getFirst() : null;
            kmerSources = akp != null ? akp.getSecond() : null;
        }

        if (!ec.getLinks().isEmpty()) {
            for (CortexLinks lm : ec.getLinks().keySet()) {
                if (lm.containsKey(curKmer)) {
                    linkStore.add(curKmer, lm.get(curKmer), false, ec.getLinks().get(lm));
                }
            }
        }

        return cv;
    }

    private Pair<ByteKmer, Set<String>> getAdjacentKmer(ByteKmer kmer, Set<CortexVertex> adjKmers, boolean goForward) {
        Pair<String, Set<String>> choicePair = linkStore.getNextJunctionChoice();
        String choice = choicePair.getFirst();
        Set<String> sources = choicePair.getSecond();

        if (choice != null) {
            Set<ByteKmer> adjKmersStr = new HashSet<>();
            for (CortexVertex cv : adjKmers) {
                adjKmersStr.add(cv.getBk());
            }

            byte[] bAdjKmer = new byte[kmer.length()];
            if (goForward) {
                System.arraycopy(kmer.getKmer(), 1, bAdjKmer, 0, kmer.length() - 1);
                bAdjKmer[bAdjKmer.length - 1] = choice.getBytes()[0];
            } else {
                bAdjKmer[0] = choice.getBytes()[0];
                System.arraycopy(kmer.getKmer(), 0, bAdjKmer, 1, kmer.length() - 1);
            }

            ByteKmer adjKmer = new ByteKmer(bAdjKmer);
            //String adjKmer = goForward ? kmer.substring(1, kmer.length()) + choice : choice + kmer.substring(0, kmer.length() - 1);

            if (adjKmersStr.contains(adjKmer)) {
                return new Pair<>(adjKmer, sources);
            }
        }

        return null;
    }

    public void setCursor(String sk, boolean goForward) {
        if (sk != null) {
            curKmer = new ByteKmer(sk.getBytes());

            Set<CortexVertex> prevKmers = getPrevVertices(curKmer);
            prevKmer = (prevKmers.size() == 1) ? prevKmers.iterator().next().getBk() : null;

            Set<CortexVertex> nextKmers = getNextVertices(curKmer);
            nextKmer = (nextKmers.size() == 1) ? nextKmers.iterator().next().getBk() : null;

            linkStore = new LinkStore();
            seen = new HashSet<>();

            if (!ec.getLinks().isEmpty()) {
                for (CortexLinks lm : ec.getLinks().keySet()) {
                    CortexKmer ck = new CortexKmer(curKmer.getKmer());
                    if (lm.containsKey(ck)) {
                        linkStore.add(curKmer, lm.get(ck), goForward, ec.getLinks().get(lm));
                    }
                }
            }
        }
    }

    //public String getCursor() { return curKmer; }

    @Override public boolean hasNext() { return nextKmer != null; }

    @Override public boolean hasPrevious() { return prevKmer != null; }

    @Override public int nextIndex() { return 0; }

    @Override public int previousIndex() { return 0; }

    @Override public void remove() { throw new UnsupportedOperationException("Cannot remove elements from CortexGraph iterator"); }

    @Override public void set(CortexVertex cortexVertex) { throw new UnsupportedOperationException("Cannot set elements in CortexGraph iterator"); }

    @Override public void add(CortexVertex cortexVertex) { throw new UnsupportedOperationException("Cannot add elements to CortexGraph iterator"); }

    /*
    private DirectedWeightedPseudograph<CortexVertex, CortexEdge> addDisplayColors(DirectedGraph<CortexVertex, CortexEdge> g) {
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> m = new DirectedWeightedPseudograph<>(CortexEdge.class);

        Set<Integer> displayColors = new HashSet<>(ec.getDisplayColors());
        if (displayColors.isEmpty()) {
            Graphs.addGraph(m, g);
        } else {
            Map<String, Map<Integer, Set<ByteKmer>>> pkscache = new HashMap<>();
            Map<String, Map<Integer, Set<ByteKmer>>> nkscache = new HashMap<>();

            for (int c : displayColors) {
                DirectedGraph<CortexVertex, CortexEdge> g2 = new DefaultDirectedGraph<>(CortexEdge.class);

                for (CortexVertex v : g.vertexSet()) {
                    Map<Integer, Set<ByteKmer>> pks = pkscache.containsKey(v.getBk()) ? pkscache.get(v.getBk()) : getAllPrevKmers(v.getBk());
                    Map<Integer, Set<ByteKmer>> nks = nkscache.containsKey(v.getBk()) ? nkscache.get(v.getBk()) : getAllNextKmers(v.getBk());

                    pkscache.put(v.getSk(), pks);
                    nkscache.put(v.getSk(), nks);

                    g2.addVertex(v);

                    for (String pk : pks.get(c)) {
                        CortexVertex pv = new CortexVertex(pk, ec.getGraph().findRecord(new CortexKmer(pk)));

                        g2.addVertex(pv);
                        g2.addEdge(pv, v, new CortexEdge(c, 1.0));
                    }

                    for (String nk : nks.get(c)) {
                        CortexVertex nv = new CortexVertex(nk, ec.getGraph().findRecord(new CortexKmer(nk)));

                        g2.addVertex(nv);
                        g2.addEdge(v, nv, new CortexEdge(c, 1.0));
                    }
                }

                Graphs.addGraph(m, g2);
            }
        }

        return m;
    }
    */
}
