package uk.ac.ox.well.indiana.commands.gg;

import com.sun.net.httpserver.Headers;
import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.json.JSONArray;
import org.json.JSONObject;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerIndex;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.*;
import java.net.InetSocketAddress;
import java.util.*;

public class VisualizeGraph extends Module {
    @Argument(fullName = "graphClean", shortName = "c", doc = "Graph (clean)")
    public CortexGraph CLEAN;

    @Argument(fullName = "graphDirty", shortName = "d", doc = "Graph (dirty)", required = false)
    public CortexGraph DIRTY;

    @Argument(fullName = "novelGraph", shortName = "n", doc = "Graph of novel kmers")
    public CortexGraph NOVEL;

    @Argument(fullName = "ref1", shortName = "r1", doc = "Fasta file for first parent")
    public File REF1;

    @Argument(fullName = "ref2", shortName = "r2", doc = "Fasta file for second parent")
    public File REF2;

    @Argument(fullName = "bed", shortName = "b", doc = "Bed file describing variants", required = false)
    public File BED;

    @Argument(fullName = "novelKmerMap", shortName = "m", doc = "Novel kmer map", required = false)
    public File NOVEL_KMER_MAP;

    @Argument(fullName="skipToKmer", shortName="s", doc="Skip processing to given kmer", required=false)
    public String KMER;

    @Argument(fullName="port", shortName="p", doc="Port")
    public Integer PORT = 9000;

    @Output
    public PrintStream out;

    private abstract class BaseHandler implements HttpHandler {
        public Map<String, String> query(String query) {
            Map<String, String> result = new HashMap<String, String>();

            query = query.replaceAll("&&", "<and>");

            for (String param : query.split("&")) {
                String pair[] = param.split("=", 2);
                if (pair.length > 1) {
                    result.put(pair[0], pair[1].replaceAll("<and>", "&&"));
                } else {
                    result.put(pair[0], "");
                }
            }

            return result;
        }

        public void write(HttpExchange httpExchange, String response) throws IOException {
            write(httpExchange, 200, response);
        }

        public void write(HttpExchange httpExchange, int code, String response) throws IOException {
            httpExchange.sendResponseHeaders(code, response.length());

            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();
        }
    }

    private class PageHandler extends BaseHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            String requestUri = httpExchange.getRequestURI().toString();
            if (requestUri.isEmpty() || requestUri.equals("/")) {
                requestUri = "index.html";
            }

            String page = "/html/" + requestUri;
            File f = new File("./public" + page);

            int code = 200;
            String response;

            if (f.exists()) {
                InputStream is;

                if (log.isDebugEnabled()) {
                    is = new FileInputStream(f);
                } else {
                    is = this.getClass().getResourceAsStream(page);
                }

                BufferedReader in = new BufferedReader(new InputStreamReader(is));
                StringBuilder responseData = new StringBuilder();

                String line;
                while ((line = in.readLine()) != null) {
                    responseData.append(line).append("\n");
                }

                Headers h = httpExchange.getResponseHeaders();
                if      (page.endsWith(".js"))   { h.add("Content-Type", "text/javascript"); }
                else if (page.endsWith(".css"))  { h.add("Content-Type", "text/css");        }
                else if (page.endsWith(".html")) { h.add("Content-Type", "text/html");       }

                response = responseData.toString();
            } else {
                code = 404;
                response = "Page not found.";
            }

            write(httpExchange, code, response);
        }
    }

    private class Graph extends BaseHandler {
        private Map<CortexKmer, Boolean> novelKmers;
        private List<CortexKmer> novelKmersList;
        private Random rng = new Random(0);
        private Map<CortexKmer, VariantInfo> vis;
        private KmerLookup kl1;
        private KmerLookup kl2;

        public Graph() {
            novelKmers = new HashMap<CortexKmer, Boolean>();
            novelKmersList = new ArrayList<CortexKmer>();

            for (CortexRecord cr : NOVEL) {
                novelKmers.put(new CortexKmer(cr.getKmerAsString()), true);
                novelKmersList.add(cr.getCortexKmer());
            }

            vis = GenotypeGraphUtils.loadNovelKmerMap(NOVEL_KMER_MAP, BED);

            kl1 = new KmerLookup(REF1);
            kl2 = new KmerLookup(REF2);
        }

        private String getRandomNovelKmer() {
            int index = rng.nextInt(novelKmersList.size());
            return novelKmersList.get(index).getKmerAsString();
        }

        private String getAppropriateKmer(String query) {
            String kmer = query;

            if (query.contains("random")) {
                boolean isOfRequestedType = false;
                String[] pieces = query.split("\\s+");
                String type = pieces.length > 1 ? pieces[1].toUpperCase() : "ANY";

                do {
                    kmer = getRandomNovelKmer();

                    CortexKmer ck = new CortexKmer(kmer);
                    VariantInfo vi = vis.containsKey(ck) ? vis.get(ck) : null;

                    if (type.equals("ANY") || (vi != null && vi.denovo.equals(type))) {
                        isOfRequestedType = true;
                    }
                } while (!isOfRequestedType);
            }

            return kmer;
        }

        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = query(httpExchange.getRequestURI().getQuery());

            /*
            String kmer = query.get("kmer");
            if (kmer.contains("random")) {
                int index = rng.nextInt(novelKmersList.size());
                kmer = novelKmersList.get(index).getKmerAsString();
            }
            */

            String kmer = getAppropriateKmer(query.get("kmer"));

            CortexKmer ck = new CortexKmer(kmer);
            VariantInfo vi = vis.containsKey(ck) ? vis.get(ck) : null;

            boolean simplify = query.get("simplify").equals("true");

            log.info("");
            log.info("Request: {}", httpExchange.getRequestURI());
            log.info("  kmer: {}", kmer);
            log.info("  known: {}", vi);

            String stretch = CortexUtils.getSeededStretch(CLEAN, DIRTY, kmer, 0, true);

            log.info("  stretch: {} bp", stretch.length());

            DirectedGraph<AnnotatedVertex, AnnotatedEdge> a = GenotypeGraphUtils.loadLocalSubgraph(stretch, CLEAN, DIRTY, novelKmers);

            //
            for (AnnotatedVertex av : a.vertexSet()) {
                if (!av.isNovel()) {
                    log.info("  aligning: {} {}", av, kl1.findKmer(av.getKmer()));

                    av.setMaternalLocations(kl1.findKmer(av.getKmer()));

                    log.info("  aligning: {} {}", av, kl2.findKmer(av.getKmer()));

                    av.setPaternalLocations(kl2.findKmer(av.getKmer()));
                }
            }
            //

            log.info("    subgraph  : {} vertices, {} edges", a.vertexSet().size(), a.edgeSet().size());

            GraphicalVariantContext gvc = new GraphicalVariantContext()
                    .attribute(0, "stretch", stretch)
                    .attribute(0, "stretchLength", stretch.length())
                    .attribute(0, "novelKmersTotal", novelKmers.size());

            // Extract parental stretches
            //
            PathInfo p1 = GenotypeGraphUtils.computeBestMinWeightPath(CLEAN, DIRTY, a, 1, stretch, novelKmers);
            PathInfo p2 = GenotypeGraphUtils.computeBestMinWeightPath(CLEAN, DIRTY, a, 2, stretch, novelKmers);

            gvc.add(GenotypeGraphUtils.callVariant(CLEAN, DIRTY, p1, 1, stretch, novelKmers, kl1));
            gvc.add(GenotypeGraphUtils.callVariant(CLEAN, DIRTY, p2, 2, stretch, novelKmers, kl2));

            // Finalize into a single call
            GenotypeGraphUtils.chooseVariant(gvc);
            //

            int numVertices = a.vertexSet().size();

            if (simplify) {
                a = GenotypeGraphUtils.simplifyGraph(a, false);

                log.info("    simplified: {} vertices, {} edges", a.vertexSet().size(), a.edgeSet().size());
            }

            int numVerticesSimplified = a.vertexSet().size();

            JSONArray va = new JSONArray();
            JSONArray ea = new JSONArray();

            Map<String, Integer> indices = new HashMap<String, Integer>();
            int index = 0;

            for (AnnotatedVertex v : a.vertexSet()) {
                Map<String, Object> m = new HashMap<String, Object>();
                m.put("kmer", v.getKmer());
                m.put("isNovel", v.isNovel());
                m.put("isPredecessor", v.flagIsSet("predecessor"));
                m.put("isSuccessor", v.flagIsSet("successor"));
                m.put("branchRejected", v.flagIsSet("branchRejected"));
                /*
                m.put("isStart", v.getKmer().equals(p1.start) || v.getKmer().equals(p2.start));
                m.put("isStop", v.getKmer().equals(p1.stop) || v.getKmer().equals(p2.stop));
                m.put("maternalLocations", v.getMaternalLocations());
                m.put("paternalLocations", v.getPaternalLocations());
                */

                va.put(m);

                indices.put(v.getKmer(), index);
                index++;
            }

            for (AnnotatedEdge e : a.edgeSet()) {
                String as = a.getEdgeSource(e).getKmer();
                String at = a.getEdgeTarget(e).getKmer();

                int is = indices.get(as);
                int it = indices.get(at);

                for (int c = 0; c < 3; c++) {
                    if (e.isPresent(c)) {
                        Map<String, Object> m = new HashMap<String, Object>();
                        m.put("source", is);
                        m.put("target", it);
                        m.put("sample", c);
                        m.put("highlight", false);

                        //
                        if (c == 0 && gvc.getAttributeAsString(0, "childStretch").contains(as) && gvc.getAttributeAsString(0, "childStretch").contains(at)) {
                            m.put("highlight", true);
                        }

                        if (c != 0 && (c == gvc.getAttributeAsInt(0, "haplotypeBackground") || gvc.getAttributeAsInt(0, "haplotypeBackground") == 0) && gvc.getAttributeAsString(0, "parentalStretch").contains(as) && gvc.getAttributeAsString(0, "parentalStretch").contains(at)) {
                            m.put("highlight", true);
                        }
                        //

                        ea.put(m);
                    }
                }
            }

            JSONObject jo = new JSONObject();
            jo.put("kmer", kmer);
            jo.put("stretch", stretch);
            jo.put("numVertices", numVertices);
            jo.put("numVerticesSimplified", numVerticesSimplified);
            jo.put("nodes", va);
            jo.put("links", ea);
            jo.put("knownVariant", vi == null ? "unknown" : vi.variantId);
            jo.put("knownType", vi == null ? "unknown" : vi.denovo);
            jo.put("knownRef", vi == null ? "unknown" : vi.ref);
            jo.put("knownAlt", vi == null ? "unknown" : vi.alt);
            jo.put("knownLeft", vi == null ? "unknown" : vi.leftFlank);
            jo.put("knownRight", vi == null ? "unknown" : vi.rightFlank);

            write(httpExchange, jo.toString());
        }
    }

    private class CrRecord extends BaseHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = query(httpExchange.getRequestURI().getQuery());

            String seq = query.get("seq");
            String sk = seq.substring(0, CLEAN.getKmerSize());
            CortexKmer ck = new CortexKmer(sk);

            boolean clean = true;
            CortexRecord cr = CLEAN.findRecord(ck);
            if (cr == null && DIRTY != null) {
                cr = DIRTY.findRecord(ck);
                clean = false;
            }

            String recordText = cr == null ? "not found" : GenotypeGraphUtils.recordToString(sk, cr) + " (" + (clean ? "clean" : "dirty") + ")";

            JSONObject jo = new JSONObject();
            jo.put("kmer", sk);
            jo.put("recordText", recordText);

            write(httpExchange, jo.toString());
        }
    }

    private class Explore extends BaseHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = query(httpExchange.getRequestURI().getQuery());

            log.info("");
            log.info("Request: {}", httpExchange.getRequestURI());
            log.info("  kmer: {}", query.get("kmer"));

            String sk = query.get("kmer");

            DirectedGraph<AnnotatedVertex, AnnotatedEdge> a = new DefaultDirectedGraph<AnnotatedVertex, AnnotatedEdge>(AnnotatedEdge.class);
            Graphs.addGraph(a, CortexUtils.dfs(CLEAN, DIRTY, sk, 0, null, ExplorationStopper.class));
            Graphs.addGraph(a, CortexUtils.dfs(CLEAN, DIRTY, sk, 1, null, ExplorationStopper.class));
            Graphs.addGraph(a, CortexUtils.dfs(CLEAN, DIRTY, sk, 2, null, ExplorationStopper.class));

            JSONArray va = new JSONArray();
            JSONArray ea = new JSONArray();

            Map<String, Integer> indices = new HashMap<String, Integer>();
            int index = Integer.valueOf(query.get("sindex"));

            for (AnnotatedVertex v : a.vertexSet()) {
                Map<String, Object> m = new HashMap<String, Object>();
                m.put("kmer", v.getKmer());
                m.put("isNovel", v.isNovel());
                m.put("isPredecessor", false);
                m.put("isSuccessor", false);
                m.put("isStart", false);
                m.put("isStop", false);
                m.put("expandedFrom", sk);

                va.put(m);

                indices.put(v.getKmer(), index);
                index++;
            }

            for (AnnotatedEdge e : a.edgeSet()) {
                String as = a.getEdgeSource(e).getKmer();
                String at = a.getEdgeTarget(e).getKmer();

                int is = indices.get(as);
                int it = indices.get(at);

                for (int c = 0; c < 3; c++) {
                    if (e.isPresent(c)) {
                        Map<String, Object> m = new HashMap<String, Object>();
                        m.put("source", is);
                        m.put("target", it);
                        m.put("sample", c);
                        m.put("highlight", false);
                        m.put("expandedFrom", sk);

                        ea.put(m);
                    }
                }
            }

            JSONObject jo = new JSONObject();
            jo.put("nodes", va);
            jo.put("links", ea);

            write(httpExchange, jo.toString());
        }
    }

    @Override
    public void execute() {
        log.info("Starting server...");

        try {
            HttpServer server = HttpServer.create(new InetSocketAddress(PORT), 0);
            server.createContext("/", new PageHandler());
            server.createContext("/graph", new Graph());
            server.createContext("/crrecord", new CrRecord());
            server.createContext("/explore", new Explore());
            server.setExecutor(null);
            server.start();
        } catch (IOException e) {
            throw new IndianaException("Unable to start server", e);
        }

        log.info("  listening on port {}", PORT);
    }
}
