package uk.ac.ox.well.cortexjdk.utils.visualizer.handlers;

import com.google.api.client.http.HttpStatusCodes;
import com.sun.net.httpserver.HttpExchange;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.json.JSONObject;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;

/**
 * Created by kiran on 29/05/2017.
 */
public class SubGraphListener extends BaseHandler {
    private Map<String, JSONObject> preLoadedGraphs = new TreeMap<>();

    public SubGraphListener() {
        DirectedGraph<CortexVertex, CortexEdge> g = new DefaultDirectedGraph<>(CortexEdge.class);

        String template = new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(100));
        StringBuilder sb1 = new StringBuilder(template);
        sb1.setCharAt(50, 'A');
        StringBuilder sb2 = new StringBuilder(template);
        sb2.setCharAt(50, 'C');

        String seq1 = sb1.toString();
        String seq2 = sb2.toString();

        String[] seqs = new String[] { seq1, seq2 };

        int kmerSize = 47;

        for (int q = 0; q < seqs.length; q++) {
            String seq = seqs[q];

            CortexVertex lv = null;

            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                String kmer = seq.substring(i, i + kmerSize);

                CortexVertex cv = new CortexVertex(new CortexByteKmer(kmer.getBytes()), null);
                g.addVertex(cv);

                if (lv != null) {
                    g.addEdge(lv, cv, new CortexEdge(q, 1.0));
                }

                lv = cv;
            }
        }

        JSONObject jo = new JSONObject();

        Set<Map<String, Object>> vs = new HashSet<>();
        Set<Map<String, Object>> es = new HashSet<>();

        for (CortexEdge e : g.edgeSet()) {
            CortexVertex v0 = g.getEdgeSource(e);
            CortexVertex v1 = g.getEdgeTarget(e);

            Map<String, Object> em = new HashMap<>();
            em.put("source", v0.getSk());
            em.put("target", v1.getSk());
            em.put("color",  e.getColor());
            es.add(em);

            Map<String, Object> va = new HashMap<>();
            va.put("id", v0.getSk());
            vs.add(va);

            Map<String, Object> vb = new HashMap<>();
            vb.put("id", v1.getSk());
            vs.add(vb);
        }

        jo.put("vertices", vs);
        jo.put("edges", es);

        preLoadedGraphs.put("test1", jo);
    }

    @Override
    public Pair<Integer, String> respond(HttpExchange httpExchange) throws IOException {
        if (httpExchange.getRequestMethod().equals("POST")) {
            StringBuilder sb = new StringBuilder();

            new BufferedReader(new InputStreamReader(httpExchange.getRequestBody()))
                    .lines()
                    .forEach((String s) -> sb.append(s).append("\n"));

            JSONObject jo = new JSONObject(sb.toString());
            preLoadedGraphs.put(jo.getString("name"), jo);
        } else if (httpExchange.getRequestMethod().equals("GET") && httpExchange.getRequestURI().getQuery() != null) {
            Map<String, String> query = query(httpExchange);

            if (query.containsKey("name")) {
                return new Pair<>(HttpStatusCodes.STATUS_CODE_OK, preLoadedGraphs.get(query.get("name")).toString());
            }
        }

        JSONObject graphNames = new JSONObject();
        graphNames.put("names", preLoadedGraphs.keySet());

        return new Pair<>(HttpStatusCodes.STATUS_CODE_OK, graphNames.toString());
    }
}
