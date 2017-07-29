package uk.ac.ox.well.indiana.utils.traversal;

import org.jetbrains.annotations.Nullable;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;
import uk.ac.ox.well.indiana.utils.stoppingconditions.TraversalStopper;

import java.util.*;
import java.util.function.Consumer;

import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.REVERSE;

public class TraversalEngine implements ListIterator<CortexVertex> {
    private TraversalEngineConfiguration ec;

    public TraversalEngine(TraversalEngineConfiguration ec) { this.ec = ec; }

    public TraversalEngineConfiguration getConfiguration() { return ec; }

    private String curKmer = null;
    private String prevKmer;
    private String nextKmer;

    private LinkStore linkStore;
    private Set<String> seen;

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

        if (dfs != null) {
            return addDisplayColors(dfs);
        }

        return null;
    }

    public String walk(String seed) {
        StringBuilder sb = new StringBuilder();
        sb.append(seed);

        setCursor(seed, true);
        while (hasNext()) {
            CortexVertex cv = next();
            String sk = cv.getSk();

            sb.append(sk.substring(sk.length() - 1, sk.length()));
        }

        setCursor(seed, false);
        while (hasPrevious()) {
            CortexVertex cv = previous();
            String sk = cv.getSk();

            sb.insert(0, sk.substring(0, 1));
        }

        return sb.toString();
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

        CortexVertex cv = new CortexVertex(sk, cr);

        Set<CortexVertex> avs;

        TraversalStopper<CortexVertex, CortexEdge> stoppingRule = instantiateStopper(ec.getStoppingRule());

        do {
            Set<CortexVertex> pvs = getPrevVertices(cv.getSk());
            Set<CortexVertex> nvs = getNextVertices(cv.getSk());
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

    public Set<CortexVertex> getPrevVertices(String sk) {
        Set<CortexVertex> prevVertices = new HashSet<>();

        Map<Integer, Set<String>> prevKmers = getAllPrevKmers(sk);

        if (prevKmers.get(ec.getTraversalColor()).size() > 0) {
            for (String prevKmer : prevKmers.get(ec.getTraversalColor())) {
                prevVertices.add(new CortexVertex(prevKmer, ec.getGraph().findRecord(new CortexKmer(prevKmer))));
            }
        } else {
            Map<String, Set<Integer>> inKmerMap = new HashMap<>();

            for (int c : ec.getRecruitmentColors()) {
                Set<String> inKmers = prevKmers.get(c);

                for (String inKmer : inKmers) {
                    if (!inKmerMap.containsKey(inKmer)) {
                        inKmerMap.put(inKmer, new HashSet<>());
                    }

                    inKmerMap.get(inKmer).add(c);
                }
            }

            for (String prevKmer : inKmerMap.keySet()) {
                prevVertices.add(new CortexVertex(prevKmer, ec.getGraph().findRecord(new CortexKmer(prevKmer))));
            }
        }

        return prevVertices;
    }

    public int getInDegree(String sk) { return getPrevVertices(sk).size(); }

    public int getOutDegree(String sk) { return getNextVertices(sk).size(); }

    public Set<CortexVertex> getNextVertices(String sk) {
        Set<CortexVertex> nextVertices = new HashSet<>();

        Map<Integer, Set<String>> nextKmers = getAllNextKmers(sk);
        if (nextKmers.get(ec.getTraversalColor()).size() > 0) {
            for (String nextKmer : nextKmers.get(ec.getTraversalColor())) {
                nextVertices.add(new CortexVertex(nextKmer, ec.getGraph().findRecord(new CortexKmer(nextKmer))));
            }
        } else {
            Map<String, Set<Integer>> outKmerMap = new HashMap<>();

            for (int c : ec.getRecruitmentColors()) {
                Set<String> outKmers = nextKmers.get(c);

                for (String outKmer : outKmers) {
                    if (!outKmerMap.containsKey(outKmer)) {
                        outKmerMap.put(outKmer, new HashSet<>());
                    }

                    outKmerMap.get(outKmer).add(c);
                }
            }

            for (String nextKmer : outKmerMap.keySet()) {
                nextVertices.add(new CortexVertex(nextKmer, ec.getGraph().findRecord(new CortexKmer(nextKmer))));
            }
        }

        return nextVertices;
    }

    public static Map<Integer, Set<String>> getAllPrevKmers(CortexRecord cr, boolean isFlipped) {
        Map<Integer, Set<String>> prevKmers = new HashMap<>();

        if (cr != null) {
            String sk = !isFlipped ? cr.getKmerAsString() : SequenceUtils.reverseComplement(cr.getKmerAsString());
            String suffix = sk.substring(0, sk.length() - 1);

            Map<Integer, Set<String>> inEdges = getInEdges(cr, isFlipped);

            for (int c = 0; c < cr.getNumColors(); c++) {
                Set<String> inKmers = new HashSet<>();

                for (String inEdge : inEdges.get(c)) {
                    inKmers.add(inEdge + suffix);
                }

                prevKmers.put(c, inKmers);
            }
        }

        return prevKmers;
    }

    public Map<Integer, Set<String>> getAllPrevKmers(String sk) {
        CortexKmer ck = new CortexKmer(sk);
        CortexRecord cr = ec.getGraph().findRecord(ck);

        return getAllPrevKmers(cr, ck.isFlipped());
    }

    public static Map<Integer, Set<String>> getAllNextKmers(CortexRecord cr, boolean isFlipped) {
        Map<Integer, Set<String>> nextKmers = new HashMap<>();

        if (cr != null) {
            String sk = !isFlipped ? cr.getKmerAsString() : SequenceUtils.reverseComplement(cr.getKmerAsString());
            String prefix = sk.substring(1, sk.length());

            Map<Integer, Set<String>> outEdges = getOutEdges(cr, isFlipped);

            for (int c = 0; c < cr.getNumColors(); c++) {
                Set<String> outKmers = new HashSet<>();

                for (String outEdge : outEdges.get(c)) {
                    outKmers.add(prefix + outEdge);
                }

                nextKmers.put(c, outKmers);
            }
        }

        return nextKmers;
    }

    public Map<Integer, Set<String>> getAllNextKmers(String sk) {
        CortexKmer ck = new CortexKmer(sk);
        CortexRecord cr = ec.getGraph().findRecord(ck);

        return getAllNextKmers(cr, ck.isFlipped());
    }

    private static Map<Integer, Set<String>> getInEdges(CortexRecord cr, boolean kmerIsFlipped) {
        Map<Integer, Set<String>> inEdges = new HashMap<>();

        if (!kmerIsFlipped) {
            for (int c = 0; c < cr.getNumColors(); c++) {
                inEdges.put(c, new HashSet<>(cr.getInEdgesAsStrings(c, false)));
            }
        } else {
            for (int c = 0; c < cr.getNumColors(); c++) {
                inEdges.put(c, new HashSet<>(cr.getOutEdgesAsStrings(c, true)));
            }
        }

        return inEdges;
    }

    private static Map<Integer, Set<String>> getOutEdges(CortexRecord cr, boolean kmerIsFlipped) {
        Map<Integer, Set<String>> outEdges = new HashMap<>();

        if (!kmerIsFlipped) {
            for (int c = 0; c < cr.getNumColors(); c++) {
                outEdges.put(c, new HashSet<>(cr.getOutEdgesAsStrings(c, false)));
            }
        } else {
            for (int c = 0; c < cr.getNumColors(); c++) {
                outEdges.put(c, new HashSet<>(cr.getInEdgesAsStrings(c, true)));
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

        CortexRecord cr = ec.getGraph().findRecord(new CortexKmer(nextKmer));
        CortexVertex cv = new CortexVertex(nextKmer, cr);

        prevKmer = curKmer;
        curKmer = nextKmer;

        Set<CortexVertex> nextKmers = getNextVertices(curKmer);
        nextKmer = null;

        if (nextKmers.size() == 1 && (!seen.contains(nextKmers.iterator().next().getSk()) || linkStore.isActive())) {
            nextKmer = nextKmers.iterator().next().getSk();

            seen.add(nextKmer);
        } else if (nextKmers.size() > 1) {
            nextKmer = getAdjacentKmer(curKmer, nextKmers, true);
        }

        if (!ec.getLinks().isEmpty()) {
            for (CortexLinks lm : ec.getLinks().keySet()) {
                if (lm.containsKey(curKmer)) {
                    linkStore.add(curKmer, lm.get(curKmer), true);
                }
            }
        }

        return cv;
    }

    @Override
    public CortexVertex previous() {
        if (prevKmer == null) { throw new NoSuchElementException("No single prev kmer from cursor '" + curKmer + "'"); }

        linkStore.incrementAge();

        CortexRecord cr = ec.getGraph().findRecord(new CortexKmer(prevKmer));
        CortexVertex cv = new CortexVertex(prevKmer, cr);

        nextKmer = curKmer;
        curKmer = prevKmer;

        Set<CortexVertex> prevKmers = getPrevVertices(curKmer);
        prevKmer = null;

        if (prevKmers.size() == 1 && (!seen.contains(prevKmers.iterator().next().getSk()) || linkStore.isActive())) {
            prevKmer = prevKmers.iterator().next().getSk();

            seen.add(prevKmer);
        } else if (prevKmers.size() > 1) {
            prevKmer = getAdjacentKmer(curKmer, prevKmers, false);
        }

        if (!ec.getLinks().isEmpty()) {
            for (CortexLinks lm : ec.getLinks().keySet()) {
                if (lm.containsKey(curKmer)) {
                    linkStore.add(curKmer, lm.get(curKmer), false);
                }
            }
        }

        return cv;
    }

    private String getAdjacentKmer(String kmer, Set<CortexVertex> adjKmers, boolean goForward) {
        String choice = linkStore.getNextJunctionChoice();

        if (choice != null) {
            Set<String> adjKmersStr = new HashSet<>();
            for (CortexVertex cv : adjKmers) {
                adjKmersStr.add(cv.getSk());
            }

            String adjKmer = goForward ? kmer.substring(1, kmer.length()) + choice : choice + kmer.substring(0, kmer.length() - 1);

            if (adjKmersStr.contains(adjKmer)) {
                return adjKmer;
            }
        }

        return null;
    }

    public void setCursor(String sk, boolean goForward) {
        curKmer = sk;

        if (curKmer != null) {
            Set<CortexVertex> prevKmers = getPrevVertices(curKmer);
            prevKmer = (prevKmers.size() == 1) ? prevKmers.iterator().next().getSk() : null;

            Set<CortexVertex> nextKmers = getNextVertices(curKmer);
            nextKmer = (nextKmers.size() == 1) ? nextKmers.iterator().next().getSk() : null;

            linkStore = new LinkStore();
            seen = new HashSet<>();

            if (!ec.getLinks().isEmpty()) {
                for (CortexLinks lm : ec.getLinks().keySet()) {
                    if (lm.containsKey(curKmer)) {
                        linkStore.add(curKmer, lm.get(curKmer), goForward);
                    }
                }
            }
        }
    }

    public String getCursor() { return curKmer; }

    @Override public boolean hasNext() { return nextKmer != null; }

    @Override public boolean hasPrevious() { return prevKmer != null; }

    @Override public int nextIndex() { return 0; }

    @Override public int previousIndex() { return 0; }

    @Override public void remove() { throw new UnsupportedOperationException("Cannot remove elements from CortexGraph iterator"); }

    @Override public void set(CortexVertex cortexVertex) { throw new UnsupportedOperationException("Cannot set elements in CortexGraph iterator"); }

    @Override public void add(CortexVertex cortexVertex) { throw new UnsupportedOperationException("Cannot add elements to CortexGraph iterator"); }

    private DirectedWeightedPseudograph<CortexVertex, CortexEdge> addDisplayColors(DirectedGraph<CortexVertex, CortexEdge> g) {
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> m = new DirectedWeightedPseudograph<>(CortexEdge.class);

        Set<Integer> displayColors = new HashSet<>(ec.getDisplayColors());
        if (displayColors.isEmpty()) {
            Graphs.addGraph(m, g);
        } else {
            Map<String, Map<Integer, Set<String>>> pkscache = new HashMap<>();
            Map<String, Map<Integer, Set<String>>> nkscache = new HashMap<>();

            for (int c : displayColors) {
                DirectedGraph<CortexVertex, CortexEdge> g2 = new DefaultDirectedGraph<>(CortexEdge.class);

                for (CortexVertex v : g.vertexSet()) {
                    Map<Integer, Set<String>> pks = pkscache.containsKey(v.getSk()) ? pkscache.get(v.getSk()) : getAllPrevKmers(v.getSk());
                    Map<Integer, Set<String>> nks = nkscache.containsKey(v.getSk()) ? nkscache.get(v.getSk()) : getAllNextKmers(v.getSk());

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

}
