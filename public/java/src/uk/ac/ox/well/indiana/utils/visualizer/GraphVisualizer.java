package uk.ac.ox.well.indiana.utils.visualizer;

import com.google.api.client.http.HttpStatusCodes;
import com.sun.net.httpserver.HttpServer;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.DirectedGraph;
import org.json.JSONObject;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.stoppingconditions.ExplorationStopper;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngineFactory;
import uk.ac.ox.well.indiana.utils.visualizer.handlers.PageHandler;
import uk.ac.ox.well.indiana.utils.visualizer.handlers.SubGraphHandler;
import uk.ac.ox.well.indiana.utils.visualizer.handlers.SubGraphListener;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.InetSocketAddress;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class GraphVisualizer {
    private GraphVisualizationConfiguration vcc;

    public GraphVisualizer(GraphVisualizationConfiguration vcc) {
        this.vcc = vcc;

        try {
            HttpServer server = HttpServer.create(new InetSocketAddress(vcc.getPort()), 0);
            server.createContext("/", new PageHandler(false));
            server.createContext("/graph", new SubGraphHandler(vcc.getGraph(), vcc.getRois()));
            server.createContext("/listener", new SubGraphListener());
            server.setExecutor(null);
            server.start();
        } catch (IOException e) {
            throw new IndianaException("Unable to start server", e);
        }

        vcc.getLogger().info("Started server on port {}...", vcc.getPort());
    }

    public GraphVisualizer(int port) {
        vcc = new GraphVisualizationConfiguration();
        vcc.setPort(port);
    }

    public int display(DirectedGraph<CortexVertex, CortexEdge> g) {
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

        jo.put("vertices", vs);
        jo.put("edges", es);

        try {
            URL obj = new URL("http://localhost:" + vcc.getPort() + "/listener");
            HttpURLConnection con = (HttpURLConnection) obj.openConnection();

            con.setRequestMethod("POST");
            con.setRequestProperty("User-Agent", "Mozilla/5.0");
            con.setRequestProperty("Accept-Language", "en-US,en;q=0.5");

            con.setDoOutput(true);
            DataOutputStream wr = new DataOutputStream(con.getOutputStream());
            wr.writeBytes(jo.toString());
            wr.flush();
            wr.close();

            int responseCode = con.getResponseCode();

            return responseCode;
        } catch (MalformedURLException e) {
            throw new IndianaException("Malformed URL", e);
        } catch (IOException e) {
            throw new IndianaException("IOException", e);
        }
    }
}
