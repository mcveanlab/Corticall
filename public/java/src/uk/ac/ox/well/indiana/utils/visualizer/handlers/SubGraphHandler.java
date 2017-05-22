package uk.ac.ox.well.indiana.utils.visualizer.handlers;

import com.sun.net.httpserver.HttpExchange;
import org.jgrapht.DirectedGraph;
import org.json.JSONObject;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created by kiran on 21/05/2017.
 */
public class SubGraphHandler extends BaseHandler {
    private DirectedGraph<CortexVertex, CortexEdge> g;

    public SubGraphHandler(DirectedGraph<CortexVertex, CortexEdge> g) {
        this.g = g;
    }

    @Override
    public void handle(HttpExchange httpExchange) throws IOException {
        JSONObject jo = new JSONObject();

        Set<Map<String, Object>> vs = new HashSet<>();
        Set<Map<String, Object>> es = new HashSet<>();

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

        jo.put("nodes", vs);
        jo.put("links", es);

        write(httpExchange, jo.toString());
    }
}
