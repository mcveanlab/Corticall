package uk.ac.ox.well.cortexjdk.utils.visualizer;

import com.sun.net.httpserver.HttpServer;
import htsjdk.samtools.reference.ReferenceSequence;
import org.jgrapht.DirectedGraph;
import org.json.JSONObject;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.visualizer.handlers.PageHandler;
import uk.ac.ox.well.cortexjdk.utils.visualizer.handlers.SubGraphHandler;
import uk.ac.ox.well.cortexjdk.utils.visualizer.handlers.SubGraphListener;

import java.io.DataOutputStream;
import java.io.IOException;
import java.net.HttpURLConnection;
import java.net.InetSocketAddress;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;

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
            throw new CortexJDKException("Unable to start server", e);
        }

        vcc.getLogger().info("Started server on port {}...", vcc.getPort());
    }

    public GraphVisualizer(int port) {
        vcc = new GraphVisualizationConfiguration();
        vcc.setPort(port);
    }

    public JSONObject display(DirectedGraph<CortexVertex, CortexEdge> g, String name) {
        JSONObject jo = new JSONObject();

        Set<Map<String, Object>> vs = new HashSet<>();
        Set<Map<String, Object>> es = new HashSet<>();

        Set<Integer> colors = new HashSet<>();
        for (CortexEdge e : g.edgeSet()) {
            colors.add(e.getColor());
        }

        for (CortexEdge e : g.edgeSet()) {
            CortexVertex s = g.getEdgeSource(e);
            CortexVertex t = g.getEdgeTarget(e);

            Map<String, Object> em = new HashMap<>();
            em.put("source", s.getSk());
            em.put("target", t.getSk());
            em.put("color", e.getColor());
            es.add(em);

            Map<String, Object> sm = new HashMap<>();
            sm.put("id", s.getSk());
            sm.put("cr", recordToString(s.getSk(), s.getCr(), colors));
            vs.add(sm);

            Map<String, Object> tm = new HashMap<>();
            tm.put("id", t.getSk());
            tm.put("cr", recordToString(t.getSk(), t.getCr(), colors));
            vs.add(tm);
        }

        jo.put("vertices", vs);
        jo.put("edges", es);
        jo.put("name", name);

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

            //return con.getResponseCode();
        } catch (MalformedURLException e) {
            throw new CortexJDKException("Malformed URL", e);
        } catch (IOException e) {
            throw new CortexJDKException("IOException", e);
        }

        return jo;
    }

    public int display(DirectedGraph<CortexVertex, CortexEdge> g, List<ReferenceSequence> refs, String name) {
        JSONObject jo = new JSONObject();

        Set<Map<String, Object>> vs = new HashSet<>();
        Set<Map<String, Object>> es = new HashSet<>();

        Set<Integer> colors = new HashSet<>();
        for (CortexEdge e : g.edgeSet()) {
            colors.add(e.getColor());
        }

        for (CortexEdge e : g.edgeSet()) {
            CortexVertex s = g.getEdgeSource(e);
            CortexVertex t = g.getEdgeTarget(e);

            Map<String, Object> em = new HashMap<>();
            em.put("source", s.getSk());
            em.put("target", t.getSk());
            em.put("color", e.getColor());
            es.add(em);

            Map<String, Object> sm = new HashMap<>();
            sm.put("id", s.getSk());
            sm.put("cr", recordToString(s.getSk(), s.getCr(), colors));
            vs.add(sm);

            Map<String, Object> tm = new HashMap<>();
            tm.put("id", t.getSk());
            tm.put("cr", recordToString(t.getSk(), t.getCr(), colors));
            vs.add(tm);
        }

        Map<String, String> contigs = new TreeMap<>();
        for (ReferenceSequence rseq : refs) {
            contigs.put(rseq.getName(), rseq.getBaseString());
        }

        jo.put("vertices", vs);
        jo.put("edges", es);
        jo.put("name", name);
        jo.put("refs", contigs);

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

            return con.getResponseCode();
        } catch (MalformedURLException e) {
            throw new CortexJDKException("Malformed URL", e);
        } catch (IOException e) {
            throw new CortexJDKException("IOException", e);
        }
    }

    private String recordToString(String sk, CortexRecord cr, Set<Integer> colors) {
        String kmer = cr.getKmerAsString();
        String cov = "";
        String ed = "";

        boolean fw = sk.equals(kmer);

        if (!fw) {
            kmer = SequenceUtils.reverseComplement(kmer);
        }

        for (int c = 0; c < cr.getNumColors(); c++) {
            if (colors.contains(c)) {
                cov += " " + c + ":" + cr.getCoverage(c);
            }
        }

        for (int c = 0; c < cr.getNumColors(); c++) {
            if (colors.contains(c)) {
                String edge = cr.getEdgesAsString(c);
                ed += " " + c + ":" + (fw ? edge : SequenceUtils.reverseComplement(edge));
            }
        }

        return kmer + " " + cov + " " + ed;
    }
}
