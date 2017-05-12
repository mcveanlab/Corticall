package uk.ac.ox.well.indiana.utils.traversal;

import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultDirectedGraph;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.util.*;

import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.REVERSE;

public class TraversalEngine {
    private TraversalEngineConfiguration ec;

    public TraversalEngine(TraversalEngineConfiguration ec) {
        this.ec = ec;

        if (this.ec.graph == null) { throw new IndianaException("Graph must be specified"); }
    }

    public DirectedGraph<CortexVertex, CortexEdge> dfs(String sk) {
        DirectedGraph<CortexVertex, CortexEdge> dfsr = (ec.td == BOTH || ec.td == REVERSE) ? dfs(sk, false, 0, new HashSet<>()) : null;
        DirectedGraph<CortexVertex, CortexEdge> dfsf = (ec.td == BOTH || ec.td == FORWARD) ? dfs(sk, true,  0, new HashSet<>()) : null;

        DirectedGraph<CortexVertex, CortexEdge> dfs = new DefaultDirectedGraph<>(CortexEdge.class);

        if (ec.gco == OR) {
            if (dfsr != null || dfsf != null) {
                if (dfsr != null) { Graphs.addGraph(dfs, dfsr); }
                if (dfsf != null) { Graphs.addGraph(dfs, dfsf); }

                return dfs;
            }
        } else {
            if (dfsr != null && dfsf != null) {
                Graphs.addGraph(dfs, dfsr);
                Graphs.addGraph(dfs, dfsf);

                return dfs;
            }
        }

        return null;
    }

    private DirectedGraph<CortexVertex, CortexEdge> dfs(String sk, boolean goForward, int currentTraversalDepth, Set<CortexVertex> history) {
        DirectedGraph<CortexVertex, CortexEdge> g = new DefaultDirectedGraph<>(CortexEdge.class);

        CortexRecord cr = ec.graph.findRecord(new CortexKmer(sk));
        CortexVertex cv = new CortexVertex(sk, cr, ec.traversalColor);

        Set<CortexVertex> avs;

        do {
            Set<CortexVertex> pvs = getPrevVertices(cv.getSk());
            Set<CortexVertex> nvs = getNextVertices(cv.getSk());
            avs = goForward ? nvs : pvs;

            if (ec.connectUnusedNeighbors) {
                connectVertex(g, cv, pvs, nvs);
            }

            if (ec.stoppingRule.keepGoing(cv, goForward, ec.traversalColor, ec.joiningColors, currentTraversalDepth, g.vertexSet().size(), avs.size(), ec.previousGraph) && !history.contains(cv)) {
                history.add(cv);

                if (avs.size() == 1) {
                    if (!ec.connectUnusedNeighbors) {
                        if (goForward) { connectVertex(g, cv, null, nvs); }
                        else { connectVertex(g, cv, pvs, null); }
                    }

                    cv = avs.iterator().next();
                } else if (avs.size() != 1) {
                    boolean childrenWereSuccessful = false;

                    for (CortexVertex av : avs) {
                        if (!av.getSk().equals(sk)) {
                            DirectedGraph<CortexVertex, CortexEdge> branch = dfs(sk, goForward, currentTraversalDepth + 1, history);

                            if (branch != null) {
                                if (!ec.connectUnusedNeighbors) {
                                    if (goForward) { connectVertex(g, cv, null, nvs); }
                                    else { connectVertex(g, cv, pvs, null); }
                                }

                                Graphs.addGraph(g, branch);
                                childrenWereSuccessful = true;
                            } else {
                                // mark a rejected traversal here
                            }
                        }
                    }

                    if (childrenWereSuccessful || ec.stoppingRule.hasTraversalSucceeded(cv, goForward, ec.traversalColor, ec.joiningColors, currentTraversalDepth, g.vertexSet().size(), avs.size(), ec.previousGraph)) {
                        return g;
                    } else {
                        // mark a rejected traversal here
                    }
                }
            } else if (ec.stoppingRule.traversalSucceeded()) {
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
                g.addEdge(pv, cv, new CortexEdge());
            }
        }

        if (nvs != null) {
            for (CortexVertex nv : nvs) {
                g.addVertex(nv);
                g.addEdge(cv, nv, new CortexEdge());
            }
        }
    }

    private Set<CortexVertex> getPrevVertices(String sk) {
        Set<CortexVertex> prevVertices = new HashSet<>();

        Map<Integer, Set<String>> prevKmers = getAllPrevKmers(sk);

        if (prevKmers.get(ec.traversalColor).size() > 0) {
            for (String prevKmer : prevKmers.get(ec.traversalColor)) {
                prevVertices.add(new CortexVertex(prevKmer, ec.graph.findRecord(new CortexKmer(prevKmer)), ec.traversalColor));
            }
        } else {
            Map<String, Set<Integer>> inKmerMap = new HashMap<>();

            for (int c : ec.recruitmentColors) {
                Set<String> inKmers = prevKmers.get(c);

                for (String inKmer : inKmers) {
                    if (!inKmerMap.containsKey(inKmer)) {
                        inKmerMap.put(inKmer, new HashSet<>());
                    }

                    inKmerMap.get(inKmer).add(c);
                }
            }

            for (String prevKmer : inKmerMap.keySet()) {
                prevVertices.add(new CortexVertex(prevKmer, ec.graph.findRecord(new CortexKmer(prevKmer)), inKmerMap.get(prevKmer)));
            }
        }

        return prevVertices;
    }

    private Set<CortexVertex> getNextVertices(String sk) {
        Set<CortexVertex> nextVertices = new HashSet<>();

        Map<Integer, Set<String>> nextKmers = getAllNextKmers(sk);
        if (nextKmers.get(ec.traversalColor).size() > 0) {
            for (String nextKmer : nextKmers.get(ec.traversalColor)) {
                nextVertices.add(new CortexVertex(nextKmer, ec.graph.findRecord(new CortexKmer(nextKmer)), ec.traversalColor));
            }
        } else {
            Map<String, Set<Integer>> outKmerMap = new HashMap<>();

            for (int c : ec.recruitmentColors) {
                Set<String> outKmers = nextKmers.get(c);

                for (String outKmer : outKmers) {
                    if (!outKmerMap.containsKey(outKmer)) {
                        outKmerMap.put(outKmer, new HashSet<>());
                    }

                    outKmerMap.get(outKmer).add(c);
                }
            }

            for (String nextKmer : outKmerMap.keySet()) {
                nextVertices.add(new CortexVertex(nextKmer, ec.graph.findRecord(new CortexKmer(nextKmer)), outKmerMap.get(nextKmer)));
            }
        }

        return nextVertices;
    }


    private Map<Integer, Set<String>> getAllPrevKmers(String sk) {
        CortexKmer ck = new CortexKmer(sk);
        CortexRecord cr = ec.graph.findRecord(ck);

        Map<Integer, Set<String>> prevKmers = new HashMap<>();

        if (cr != null) {
            String suffix = sk.substring(0, sk.length() - 1);

            Map<Integer, Set<String>> inEdges = getInEdges(cr, ck.isFlipped());

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

    private Map<Integer, Set<String>> getAllNextKmers(String sk) {
        CortexKmer ck = new CortexKmer(sk);
        CortexRecord cr = ec.graph.findRecord(ck);

        Map<Integer, Set<String>> nextKmers = new HashMap<>();

        if (cr != null) {
            String prefix = sk.substring(1, sk.length());

            Map<Integer, Set<String>> outEdges = getOutEdges(cr, ck.isFlipped());

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

    private Map<Integer, Set<String>> getInEdges(CortexRecord cr, boolean kmerIsFlipped) {
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

    private Map<Integer, Set<String>> getOutEdges(CortexRecord cr, boolean kmerIsFlipped) {
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
}
