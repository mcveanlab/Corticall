package uk.ac.ox.well.indiana.utils.traversal;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.stoppingconditions.TraversalStopper;

import java.util.Set;

public class TraversalEngine {
    private TraversalEngineConfiguration configuration;

    public TraversalEngine(TraversalEngineConfiguration configuration) { this.configuration = configuration; }

    public DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs(String kmer) {



        /*
        String firstKmer = new String(kmer);

        Map<Integer, Set<String>> sourceKmersAllColors = goForward ? CortexUtils.getPrevKmers(graph, dirty, kmer) : CortexUtils.getNextKmers(graph, dirty, kmer);
        Set<String> sourceKmers = sourceKmersAllColors.get(color);

        DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs = new DefaultDirectedGraph<>(AnnotatedEdge.class);
        TraversalStopper<AnnotatedVertex, AnnotatedEdge> stopper = instantiateStopper(stopperClass);

        Map<Integer, Set<String>> adjKmers;

        Set<Integer> allColors = new HashSet<>();
        allColors.add(color);
        if (parentColors != null) { allColors.addAll(parentColors); }

        do {
            AnnotatedVertex cv = new AnnotatedVertex(kmer);

            CortexRecord cr = graph.findRecord(new CortexKmer(cv.getKmer()));
            if (cr == null && dirty != null) {
                cr = dirty.findRecord(new CortexKmer(cv.getKmer()));
            }

            Map<Integer, Set<String>> prevKmers = CortexUtils.getPrevKmers(graph, dirty, cv.getKmer());
            Map<Integer, Set<String>> nextKmers = CortexUtils.getNextKmers(graph, dirty, cv.getKmer());
            adjKmers = goForward ? nextKmers : prevKmers;

            int numVerticesAdded = addVertexAndConnect(dfs, cv, prevKmers, nextKmers, allColors);

            if (stopper.keepGoing(cr, g, depth, dfs.vertexSet().size(), adjKmers.get(color).size(), color, parentColors) && !sourceKmers.contains(kmer) && !history.contains(kmer)) {
                history.add(kmer);

                if (adjKmers.get(color).size() == 1) {
                    kmer = adjKmers.get(color).iterator().next();
                } else if (adjKmers.get(color).size() != 1) {
                    boolean childrenWereSuccessful = false;

                    for (String ak : adjKmers.get(color)) {
                        if (!ak.equals(firstKmer)) {
                            DirectedGraph<AnnotatedVertex, AnnotatedEdge> branch = dfs(graph, dirty, ak, color, parentColors, g, stopperClass, depth + (CortexUtils.isNovelKmer(cr) ? 0 : 1), goForward, history);

                            if (branch != null) {
                                Graphs.addGraph(dfs, branch);
                                childrenWereSuccessful = true;
                            } else {
                                for (AnnotatedVertex av : dfs.vertexSet()) {
                                    if (av.getKmer().equals(ak)) {
                                        av.setFlag("branchRejected");
                                    }
                                }
                            }
                        }
                    }

                    if (childrenWereSuccessful || stopper.hasTraversalSucceeded(cr, g, depth, dfs.vertexSet().size(), 0, color, parentColors)) {
                        return dfs;
                    } else {
                        for (AnnotatedVertex av : dfs.vertexSet()) {
                            if (av.getKmer().equals(kmer)) {
                                av.setFlag("branchRejected");
                            }
                        }
                    }
                }
            } else if (stopper.traversalSucceeded()) {
                return dfs;
            } else {
                return null;
            }
        } while (adjKmers.get(color).size() == 1);
        */

        return null;
    }
}
