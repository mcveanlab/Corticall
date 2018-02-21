package uk.ac.ox.well.cortexjdk.utils.traversal;

import org.jgrapht.Graph;
import org.jgrapht.GraphPath;
import org.jgrapht.alg.shortestpath.KShortestPaths;
import org.jgrapht.graph.DefaultDirectedGraph;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by kiran on 03/06/2017.
 */
public class PathFinder {
    private Graph<CortexVertex, CortexEdge> g;

    public PathFinder(Graph<CortexVertex, CortexEdge> graph, int color) {
        g = new DefaultDirectedGraph<>(CortexEdge.class);

        for (CortexEdge e : graph.edgeSet()) {
            if (e.getColor() == color) {
                CortexVertex s = graph.getEdgeSource(e);
                CortexVertex t = graph.getEdgeTarget(e);

                g.addVertex(s);
                g.addVertex(t);
                g.addEdge(s, t, e);
            }
        }
    }

    public GraphPath<CortexVertex, CortexEdge> getPath(CortexVertex startVertex, CortexVertex endVertex) {
        return getPath(startVertex, endVertex, null, true);
    }

    public GraphPath<CortexVertex, CortexEdge> getPath(CortexVertex startVertex, CortexVertex endVertex, CanonicalKmer constraint, boolean accept) {
        List<GraphPath<CortexVertex, CortexEdge>> pathsFiltered = getPaths(startVertex, endVertex, constraint, accept);

        return pathsFiltered.size() == 0 ? null : pathsFiltered.get(0);
    }

    public List<GraphPath<CortexVertex, CortexEdge>> getPaths(CortexVertex startVertex, CortexVertex endVertex) {
        return getPaths(startVertex, endVertex, null, true);
    }

    public List<GraphPath<CortexVertex, CortexEdge>> getPaths(CortexVertex startVertex, CortexVertex endVertex, CanonicalKmer constraint, boolean accept) {
        KShortestPaths<CortexVertex, CortexEdge> ksp = new KShortestPaths<>(g, 10, 500);

        List<GraphPath<CortexVertex, CortexEdge>> pathsFiltered = new ArrayList<>();

        if (g.containsVertex(startVertex) && g.containsVertex(endVertex)) {
            List<GraphPath<CortexVertex, CortexEdge>> pathsUnfiltered = ksp.getPaths(startVertex, endVertex);

            if (constraint == null) {
                pathsFiltered = pathsUnfiltered;
            } else {
                for (GraphPath<CortexVertex, CortexEdge> gp : pathsUnfiltered) {
                    boolean constraintFound = false;

                    for (CortexVertex cv : gp.getVertexList()) {
                        CanonicalKmer ck = new CanonicalKmer(cv.getKmerAsString());

                        if (ck.equals(constraint)) {
                            constraintFound = true;
                            break;
                        }
                    }

                    if ((constraintFound && accept) || (!constraintFound && !accept)) {
                        pathsFiltered.add(gp);
                    }
                }
            }
        }

        return pathsFiltered;
    }
}
