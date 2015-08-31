package uk.ac.ox.well.indiana.commands.gg;

import com.sun.net.httpserver.Headers;
import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;
import htsjdk.variant.variantcontext.VariantContext;
import org.jgrapht.DirectedGraph;
import org.json.JSONArray;
import org.json.JSONObject;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksMap;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.*;
import java.net.InetSocketAddress;
import java.util.*;

public class VisualizeGraph extends Module {
    @Argument(fullName = "graph", shortName = "g", doc = "Graph")
    public CortexGraph GRAPH;

    @Argument(fullName = "graphRaw", shortName = "r", doc = "Graph (raw)", required = false)
    public CortexGraph GRAPH_RAW;

    @Argument(fullName = "links", shortName = "l", doc = "Links", required=false)
    public CortexLinksMap LINKS;

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

    @Argument(fullName = "beAggressive", shortName = "a", doc = "Be aggressive in extending novel stretches")
    public Boolean AGGRESSIVE = false;

    @Argument(fullName = "maxJunctionsAllowed", shortName = "j", doc = "Maximum number of junctions we'll permit ourselves to traverse")
    public Integer MAX_JUNCTIONS_ALLOWED = 10;

    @Argument(fullName="skipToKmer", shortName="s", doc="Skip processing to given kmer", required=false)
    public String KMER;

    @Argument(fullName="port", shortName="p", doc="Port")
    public Integer PORT = 9000;

    @Output(fullName = "gout", shortName = "go", doc = "Graph out")
    public File gout;

    @Output
    public PrintStream out;

    @Output(fullName = "evalOut", shortName = "eout", doc = "Eval out")
    public PrintStream eout;

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

        public Graph() {
            novelKmers = new HashMap<CortexKmer, Boolean>();
            novelKmersList = new ArrayList<CortexKmer>();

            for (CortexRecord cr : NOVEL) {
                novelKmers.put(new CortexKmer(cr.getKmerAsString()), true);
                novelKmersList.add(cr.getCortexKmer());
            }

            vis = GenotypeGraphUtils.loadNovelKmerMap(NOVEL_KMER_MAP, BED);
        }

        private DirectedGraph<AnnotatedVertex, AnnotatedEdge> fetchGraph(String stretch) {
            return GenotypeGraphUtils.dfsGraph(stretch, GRAPH, GRAPH_RAW, LINKS, AGGRESSIVE, novelKmers);
        }

        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = query(httpExchange.getRequestURI().getQuery());

            String kmer = query.get("kmer");
            if (kmer.contains("random")) {
                int index = rng.nextInt(novelKmersList.size());
                kmer = novelKmersList.get(index).getKmerAsString();
            }

            CortexKmer ck = new CortexKmer(kmer);

            VariantInfo vi = vis.containsKey(ck) ? vis.get(ck) : null;

            boolean simplify = query.get("simplify").equals("true");

            log.info("");
            log.info("Request: {}", httpExchange.getRequestURI());

            String stretch = CortexUtils.getSeededStretch(GRAPH, GRAPH_RAW, kmer, 0, true);

            log.info("  stretch: {} bp", stretch.length());

            DirectedGraph<AnnotatedVertex, AnnotatedEdge> a = fetchGraph(stretch);

            log.info("    subgraph  : {} vertices, {} edges", a.vertexSet().size(), a.edgeSet().size());

            int numVertices = a.vertexSet().size();

            if (simplify) {
                a = GenotypeGraphUtils.simplifyGraph(a, false);

                log.info("    simplified: {} vertices, {} edges", a.vertexSet().size(), a.edgeSet().size());
            }

            int numVerticesSimplified = a.vertexSet().size();

            String seq = vi.leftFlank + vi.ref + vi.rightFlank;

            for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                String fw = seq.substring(i, i + GRAPH.getKmerSize());
                String rc = SequenceUtils.reverseComplement(fw);

                AnnotatedVertex af = new AnnotatedVertex(fw);
                AnnotatedVertex ar = new AnnotatedVertex(rc);

                CortexKmer ca = new CortexKmer(fw);
                CortexRecord cr = GRAPH.findRecord(ca);

                if (cr != null) {
                    log.info("  ar: clean {} {} {}", a.vertexSet().contains(ar), rc, cr);
                } else {
                    cr = GRAPH_RAW.findRecord(ca);

                    if (cr != null) {
                        log.info("  ar: dirty {} {} {}", a.vertexSet().contains(ar), rc, cr);
                    } else {
                        log.info("  ar: unknown {} {}", a.vertexSet().contains(ar), rc);
                    }
                }
            }

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
        private String recordToString(String sk, CortexRecord cr) {
            String kmer = cr.getKmerAsString();
            String cov = "";
            String ed = "";

            boolean fw = sk.equals(kmer);

            if (!fw) {
                kmer = SequenceUtils.reverseComplement(kmer);
            }

            for (int coverage : cr.getCoverages()) {
                cov += " " + coverage;
            }

            for (String edge : cr.getEdgeAsStrings()) {
                ed += " " + (fw ? edge : SequenceUtils.reverseComplement(edge));
            }

            return kmer + " " + cov + " " + ed;
        }

        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = query(httpExchange.getRequestURI().getQuery());

            String seq = query.get("seq");
            String sk = seq.substring(0, GRAPH.getKmerSize());
            CortexKmer ck = new CortexKmer(sk);

            boolean clean = true;
            CortexRecord cr = GRAPH.findRecord(ck);
            if (cr == null && GRAPH_RAW != null) {
                cr = GRAPH_RAW.findRecord(ck);
                clean = false;
            }

            String recordText = cr == null ? "not found" : recordToString(sk, cr) + " (" + (clean ? "clean" : "dirty") + ")";

            JSONObject jo = new JSONObject();
            jo.put("kmer", sk);
            jo.put("recordText", recordText);

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
            server.setExecutor(null);
            server.start();
        } catch (IOException e) {
            throw new IndianaException("Unable to start server", e);
        }

        log.info("  listening on port {}", PORT);
    }
}
