package uk.ac.ox.well.cortexjdk.utils.visualizer.handlers;

import com.google.api.client.http.HttpStatusCodes;
import com.sun.net.httpserver.HttpExchange;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.DirectedGraph;
import org.json.JSONObject;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.stoppingconditions.ExplorationStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created by kiran on 21/05/2017.
 */
public class SubGraphHandler extends BaseHandler {
    private DirectedGraph<CortexVertex, CortexEdge> g;
    private CortexGraph graph;
    private CortexGraph rois;

    public SubGraphHandler(DirectedGraph<CortexVertex, CortexEdge> g) {
        this.g = g;
    }

    public SubGraphHandler(CortexGraph graph, CortexGraph rois) {
        this.graph = graph;
        this.rois = rois;
    }

    public Pair<Integer, String> respond(HttpExchange httpExchange) {
        JSONObject jo = new JSONObject();

        Set<Map<String, Object>> vs = new HashSet<>();
        Set<Map<String, Object>> es = new HashSet<>();

        Map<String, String> qm = query(httpExchange);
        String kmer = qm.getOrDefault("kmer", null);

        if (kmer != null) {
            TraversalEngine e = new TraversalEngineFactory()
                    .graph(graph)
                    .rois(rois)
                    .traversalColor(0)
                    .stoppingRule(ExplorationStopper.class)
                    .make();
        }

        for (CortexVertex v : g.vertexSet()) {
            Map<String, Object> vm = new HashMap<>();

            vm.put("id", v.getSk());
            vs.add(vm);
        }

        for (CortexEdge e : g.edgeSet()) {
            Map<String, Object> em = new HashMap<>();

            em.put("source", g.getEdgeSource(e).getSk());
            em.put("target", g.getEdgeTarget(e).getSk());
            em.put("color", e.getColor());
            es.add(em);
        }

        jo.put("vertices", vs);
        jo.put("edges", es);

        return new Pair<>(HttpStatusCodes.STATUS_CODE_OK, jo.toString());
    }
}
