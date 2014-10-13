package uk.ac.ox.well.indiana.commands.cortex;

import com.google.common.base.Joiner;
import com.sun.net.httpserver.Headers;
import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.json.JSONArray;
import org.json.JSONObject;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexJunctionsRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksRecord;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.*;
import java.net.InetSocketAddress;
import java.util.*;

@Description(text="Starts the server for visualizing assembly data")
public class server extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs (FASTA)")
    public File CONTIGS;

    @Argument(fullName="graph", shortName="g", doc="Cortex graph file")
    public TreeMap<String, File> GRAPHS;

    @Argument(fullName="link", shortName="l", doc="Link information")
    public TreeMap<String, File> LINKS;

    @Argument(fullName="port", shortName="p", doc="Port")
    public Integer PORT = 9000;

    private Map<String, String> contigs;
    private Map<String, Map<CortexKmer, CortexLinksRecord>> links;

    private class PageHandler implements HttpHandler {
        private String page;

        public PageHandler(String page) {
            this.page = page;
        }

        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            InputStream is;

            if (log.isDebugEnabled()) {
                is = new FileInputStream(new File("./public" + page));

                log.debug("Reading '{}' from disk", page);
            } else {
                is = this.getClass().getResourceAsStream(page);
            }

            BufferedReader in = new BufferedReader(new InputStreamReader(is));
            StringBuilder responseData = new StringBuilder();

            String line;
            while((line = in.readLine()) != null) {
                responseData.append(line).append("\n");
            }

            String response = responseData.toString();
            response = response.replaceAll("\\$NUM_CONTIGS", String.valueOf(contigs.size()));

            if (page.endsWith(".js")) {
                Headers h = httpExchange.getResponseHeaders();
                h.add("Content-Type", "text/javascript");
            }

            httpExchange.sendResponseHeaders(200, response.length());

            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("Fetch static page : {}; response length: {}", httpExchange.getRequestURI(), response.length());
        }
    }

    private class ContigsHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            StringBuilder rb = new StringBuilder();
            rb.append("\"contigName\",\"baseLength\"\n");
            for (String name : contigs.keySet()) {
                rb.append("\"").append(name).append("\",\"").append(contigs.get(name).length()).append("\"\n");
            }

            String response = rb.toString();

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("GET contigs list  : {}; response length: {}", httpExchange.getRequestURI(), response.length());
        }
    }

    private Map<String, String> queryToMap(String query){
        Map<String, String> result = new HashMap<String, String>();
        for (String param : query.split("&")) {
            String pair[] = param.split("=");
            if (pair.length>1) {
                result.put(pair[0], pair[1]);
            }else{
                result.put(pair[0], "");
            }
        }
        return result;
    }

    private class ContigHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = queryToMap(httpExchange.getRequestURI().getQuery());

            StringBuilder rb = new StringBuilder();
            rb.append("{");
            rb.append("\"seq\": \"").append(contigs.get(query.get("contigName"))).append("\"");
            rb.append("}\n");

            String response = rb.toString();

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("GET contig seq    : {}, response length: {}", query.get("contigName"), response.length());
        }
    }

    private class EdgeHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = queryToMap(httpExchange.getRequestURI().getQuery());

            String seq = contigs.get(query.get("contigName"));
            String selectedGraph = query.get("graphName");

            Map<Integer, Set<String>> inEdges = new HashMap<Integer, Set<String>>();
            Map<Integer, Set<String>> outEdges = new HashMap<Integer, Set<String>>();

            for (String graphLabel : GRAPHS.keySet()) {
                if (selectedGraph.equals("all") || graphLabel.equals(selectedGraph)) {
                    CortexGraph cg = new CortexGraph(GRAPHS.get(graphLabel));
                    int kmerSize = cg.getKmerSize();

                    for (int i = 0; i <= seq.length() - kmerSize; i++) {
                        int iePos = i - 1;
                        int oePos = i + kmerSize;

                        String kmer = seq.substring(i, i + kmerSize);
                        CortexKmer ck = new CortexKmer(kmer);
                        CortexRecord cr = cg.findRecord(ck);

                        if (cr != null) {
                            Set<String> ie = CortexUtils.getPrevKmers(cg, kmer);
                            Set<String> oe = CortexUtils.getNextKmers(cg, kmer);

                            if (i > 0) {
                                String prevKmer = seq.substring(iePos, oePos - 1);
                                ie.remove(prevKmer);
                            }

                            if (i < seq.length() - kmerSize) {
                                String nextKmer = seq.substring(i + 1, i + 1 + kmerSize);
                                oe.remove(nextKmer);
                            }

                            if (!inEdges.containsKey(iePos)) {
                                inEdges.put(iePos, new TreeSet<String>());
                            }
                            if (!outEdges.containsKey(oePos)) {
                                outEdges.put(oePos, new TreeSet<String>());
                            }

                            for (String inKmer : ie) {
                                inEdges.get(iePos).add(inKmer.substring(0, 1));
                            }

                            for (String outKmer : oe) {
                                outEdges.get(oePos).add(outKmer.substring(outKmer.length() - 1, outKmer.length()));
                            }
                        }
                    }
                }
            }

            StringBuilder rb = new StringBuilder();

            boolean pastFirstRecord = false;

            rb.append("{");
            for (int i = 0; i < seq.length(); i++) {
                if ((inEdges.containsKey(i) && inEdges.get(i).size() > 0) || (outEdges.containsKey(i) && outEdges.get(i).size() > 0)) {
                    if (pastFirstRecord) {
                        rb.append(",\n\t\"").append(i).append("\": { ");
                    } else {
                        rb.append("\n\t\"").append(i).append("\": { ");
                        pastFirstRecord = true;
                    }

                    List<String> edgeArrays = new ArrayList<String>();

                    if (inEdges.containsKey(i) && inEdges.get(i).size() > 0) {
                        List<String> bases = new ArrayList<String>();
                        for (String inEdge : inEdges.get(i)) {
                            bases.add("\"" + inEdge + "\"");
                        }

                        edgeArrays.add("\"in\": [ " + Joiner.on(", ").join(bases) + " ]");
                    }

                    if (outEdges.containsKey(i) && outEdges.get(i).size() > 0) {
                        List<String> bases = new ArrayList<String>();
                        for (String outEdge : outEdges.get(i)) {
                            bases.add("\"" + outEdge + "\"");
                        }

                        edgeArrays.add("\"out\": [ " + Joiner.on(", ").join(bases) + " ]");
                    }

                    rb.append(Joiner.on(", ").join(edgeArrays)).append(" }");
                }
            }
            rb.append("\n}\n");

            String response = rb.toString();

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("GET in/out edges  : {}, response length: {}", query.get("contigName"), response.length());
        }
    }

    private class MaskHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = queryToMap(httpExchange.getRequestURI().getQuery());

            String seq = contigs.get(query.get("contigName"));
            String selectedGraph = query.get("graphName");

            Set<Integer> maskedPositions = new TreeSet<Integer>();

            for (String graphLabel : GRAPHS.keySet()) {
                if (selectedGraph.equals("all") || graphLabel.equals(selectedGraph)) {
                    CortexGraph cg = new CortexGraph(GRAPHS.get(graphLabel));
                    int kmerSize = cg.getKmerSize();

                    for (int i = 0; i <= seq.length() - kmerSize; i++) {
                        String thisKmer = seq.substring(i, i + kmerSize);
                        CortexKmer thisCk = new CortexKmer(thisKmer);
                        CortexRecord thisCr = cg.findRecord(thisCk);

                        if (thisCr == null) {
                            maskedPositions.add(i);
                        }
                    }
                }
            }

            JSONObject jsonResponse = new JSONObject();
            jsonResponse.put("masked", new JSONArray(maskedPositions));

            String response = jsonResponse.toString();

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("GET masked bases  : current, response length: {}", response.length());
        }
    }

    private class LinksHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = queryToMap(httpExchange.getRequestURI().getQuery());

            String seq = contigs.get(query.get("contigName"));
            String selectedGraph = query.get("graphName");
            boolean goForward = query.get("orientation").equals("forward");

            Set<Integer> linkStarts = new TreeSet<Integer>();

            for (String graphLabel : LINKS.keySet()) {
                if (selectedGraph.equals("all") || graphLabel.equals(selectedGraph)) {
                    Map<CortexKmer, CortexLinksRecord> cprs = links.get(graphLabel);
                    CortexGraph cg = new CortexGraph(GRAPHS.get(graphLabel));
                    int kmerSize = cg.getKmerSize();

                    for (int i = 0; i <= seq.length() - kmerSize; i++) {
                        String thisKmer = seq.substring(i, i + kmerSize);
                        CortexKmer thisCk = new CortexKmer(thisKmer);

                        if (cprs.containsKey(thisCk)) {
                            CortexLinksRecord clr = cprs.get(thisCk);

                            for (CortexJunctionsRecord cjr : clr.getJunctions()) {
                                if (thisCk.isFlipped()) {
                                    cjr = CortexUtils.flipJunctionsRecord(cjr);
                                }

                                if (goForward == cjr.isForward()) {
                                    if (goForward) {
                                        linkStarts.add(i);
                                    } else {
                                        linkStarts.add(i + kmerSize - 1);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            JSONObject jsonResponse = new JSONObject();
            jsonResponse.put("linkStarts", new JSONArray(linkStarts));

            String response = jsonResponse.toString();

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("GET links info    : current, response length: {}", response.length());
        }
    }

    private class LinkEntry {
        public String linkStartId;
        public String linkEndId;

        public LinkEntry(String linkStartId, String linkEndId) {
            this.linkStartId = linkStartId;
            this.linkEndId = linkEndId;
        }

        public Map<String, String> getJSONMap() {
            Map<String, String> jsonMap = new LinkedHashMap<String, String>();
            jsonMap.put("linkStartId", linkStartId);
            jsonMap.put("linkEndId", linkEndId);

            return jsonMap;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            LinkEntry linkEntry = (LinkEntry) o;

            if (linkEndId != null ? !linkEndId.equals(linkEntry.linkEndId) : linkEntry.linkEndId != null) return false;
            if (linkStartId != null ? !linkStartId.equals(linkEntry.linkStartId) : linkEntry.linkStartId != null)
                return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result = linkStartId != null ? linkStartId.hashCode() : 0;
            result = 31 * result + (linkEndId != null ? linkEndId.hashCode() : 0);
            return result;
        }
    }

    private class LinkHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = queryToMap(httpExchange.getRequestURI().getQuery());

            String seq = contigs.get(query.get("contigName"));
            String selectedGraph = query.get("graphName");
            int pos = Integer.valueOf(query.get("pos"));

            //Set<LinkEntry> linkEntries = new HashSet<LinkEntry>();
            Map<LinkEntry, Integer> linkEntries = new LinkedHashMap<LinkEntry, Integer>();

            for (String graphLabel : LINKS.keySet()) {
                if (selectedGraph.equals("all") || graphLabel.equals(selectedGraph)) {
                    Map<CortexKmer, CortexLinksRecord> clrs = links.get(graphLabel);
                    CortexGraph cg = new CortexGraph(GRAPHS.get(graphLabel));
                    int kmerSize = cg.getKmerSize();

                    String sk = seq.substring(pos, pos + kmerSize);
                    CortexKmer ck = new CortexKmer(sk);

                    if (!clrs.containsKey(ck)) {
                        sk = seq.substring(pos - kmerSize + 1, pos + 1);
                        ck = new CortexKmer(sk);
                    }

                    for (CortexJunctionsRecord cjr : clrs.get(ck).getJunctions()) {
                        if (ck.isFlipped()) {
                            cjr = CortexUtils.flipJunctionsRecord(cjr);
                        }

                        String junctions = cjr.getJunctions();

                        if (cjr.isForward()) {
                            log.info("junctions: {}", cjr);

                            int currentJunction = 0;
                            int startPos = pos;
                            String startBase = sk.substring(0, 1);

                            for (int i = pos; i <= seq.length() - kmerSize && currentJunction < junctions.length(); i++) {
                                String thisKmer = seq.substring(i, i + kmerSize);
                                Set<String> nextKmers = CortexUtils.getNextKmers(cg, thisKmer);
                                if (nextKmers.size() > 1) {
                                    String expectedNextKmer = seq.substring(i + 1, i + kmerSize) + junctions.charAt(currentJunction);

                                    if (nextKmers.contains(expectedNextKmer)) {
                                        int endPos = i + kmerSize;
                                        String endBase = junctions.substring(currentJunction, currentJunction + 1);

                                        LinkEntry le1 = new LinkEntry(
                                            String.format("%s_%d", startBase, startPos),
                                            String.format("%s_%d", expectedNextKmer.substring(expectedNextKmer.length() - 2, expectedNextKmer.length() - 1), i + kmerSize - 1));

                                        LinkEntry le2 = new LinkEntry(
                                            String.format("%s_%d", expectedNextKmer.substring(expectedNextKmer.length() - 2, expectedNextKmer.length() - 1), i + kmerSize - 1),
                                            String.format("%s_%d", endBase, endPos));

                                        if (linkEntries.containsKey(le1)) {
                                            linkEntries.put(le1, linkEntries.get(le1) + cjr.getCoverage(0));
                                        } else {
                                            linkEntries.put(le1, cjr.getCoverage(0));
                                        }

                                        if (linkEntries.containsKey(le2)) {
                                            linkEntries.put(le2, linkEntries.get(le2) + cjr.getCoverage(0));
                                        } else {
                                            linkEntries.put(le2, cjr.getCoverage(0));
                                        }

                                        startPos = endPos;
                                        startBase = endBase;

                                        currentJunction++;

                                        if (i < seq.length() - kmerSize) {
                                            String nextKmerInContig = seq.substring(i + 1, i + 1 + kmerSize);

                                            if (!expectedNextKmer.equals(nextKmerInContig)) {
                                                break;
                                            }
                                        }
                                    } else {
                                        log.info("Did not find expected next kmer: {} {}", i, expectedNextKmer);
                                    }
                                }
                            }
                        } else {
                            log.info("junctions: {}", cjr);

                            int currentJunction = 0;
                            int startPos = pos;
                            String startBase = sk.substring(sk.length() - 1);

                            for (int i = pos; i >= kmerSize && currentJunction < junctions.length(); i--) {
                                String thisKmer = seq.substring(i - kmerSize + 1, i + 1);
                                Set<String> prevKmers = CortexUtils.getPrevKmers(cg, thisKmer);
                                if (prevKmers.size() > 1) {
                                    String expectedPrevKmer = junctions.charAt(currentJunction) + seq.substring(i - kmerSize + 1, i);

                                    if (prevKmers.contains(expectedPrevKmer)) {
                                        int endPos = i - kmerSize;
                                        String endBase = junctions.substring(currentJunction, currentJunction + 1);

                                        LinkEntry le1 = new LinkEntry(
                                                String.format("%s_%d", startBase, startPos),
                                                String.format("%s_%d", expectedPrevKmer.substring(1, 2), i - kmerSize + 1));

                                        LinkEntry le2 = new LinkEntry(
                                                String.format("%s_%d", expectedPrevKmer.substring(1, 2), i - kmerSize + 1),
                                                String.format("%s_%d", endBase, endPos));

                                        if (linkEntries.containsKey(le1)) {
                                            linkEntries.put(le1, linkEntries.get(le1) + cjr.getCoverage(0));
                                        } else {
                                            linkEntries.put(le1, cjr.getCoverage(0));
                                        }

                                        if (linkEntries.containsKey(le2)) {
                                            linkEntries.put(le2, linkEntries.get(le2) + cjr.getCoverage(0));
                                        } else {
                                            linkEntries.put(le2, cjr.getCoverage(0));
                                        }

                                        startPos = endPos;
                                        startBase = endBase;

                                        currentJunction++;

                                        if (i > kmerSize) {
                                            String prevKmerInContig = seq.substring(i - kmerSize, i);

                                            if (!expectedPrevKmer.equals(prevKmerInContig)) {
                                                break;
                                            }
                                        }
                                    } else {
                                        log.info("Did not find expected prev kmer: {} {}", i, expectedPrevKmer);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            JSONObject jsonResponse = new JSONObject();

            for (LinkEntry le : linkEntries.keySet()) {
                Map<String, String> jsonMap = le.getJSONMap();
                jsonMap.put("linkCoverage", String.valueOf(linkEntries.get(le)));

                jsonResponse.append("linkEntries", jsonMap);
            }

            String response = jsonResponse.toString();

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("GET link info     : current, response length: {}", response.length());
        }
    }

    private class DataSetHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            JSONObject jsonResponse = new JSONObject();
            jsonResponse.put("contigsFile", CONTIGS.getAbsolutePath());

            JSONArray graphInfos = new JSONArray();
            for (String graphName : GRAPHS.keySet()) {
                Map<String, Object> graphInfo = new HashMap<String, Object>();
                CortexGraph cg = new CortexGraph(GRAPHS.get(graphName));

                graphInfo.put("graphName", graphName);
                graphInfo.put("graphFile", GRAPHS.get(graphName).getAbsolutePath());
                graphInfo.put("kmerSize", cg.getKmerSize());
                graphInfo.put("numRecords", cg.getNumRecords());
                graphInfo.put("sampleName", cg.getColor(0).getSampleName());

                graphInfos.put(graphInfo);
            }

            jsonResponse.put("graphs", graphInfos);

            String response = jsonResponse.toString();

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("GET dataset info  : current, response length: {}", response.length());
        }
    }

    private Map<String, String> loadContigs() {
        FastaSequenceFile contigReader = new FastaSequenceFile(CONTIGS, true);

        Map<String, String> contigs = new LinkedHashMap<String, String>();

        ReferenceSequence rseq;
        while ((rseq = contigReader.nextSequence()) != null) {
            String name = rseq.getName();

            contigs.put(name, new String(rseq.getBases()));
        }

        return contigs;
    }

    private Map<String, Map<CortexKmer, CortexLinksRecord>> loadLinks() {
        Map<String, Map<CortexKmer, CortexLinksRecord>> links = new HashMap<String, Map<CortexKmer, CortexLinksRecord>>();

        for (String graphName : LINKS.keySet()) {
            CortexLinks cp = new CortexLinks(LINKS.get(graphName));

            links.put(graphName, new HashMap<CortexKmer, CortexLinksRecord>());

            for (CortexLinksRecord cpr : cp) {
                links.get(graphName).put(cpr.getKmer(), cpr);
            }
        }

        return links;
    }


    @Override
    public void execute() {
        log.info("Loading...");

        contigs = loadContigs();
        links = loadLinks();

        log.info("  loaded {} contigs", contigs.size());
        for (String graphLabel : links.keySet()) {
            log.info("  loaded {} kmers with links from {}", links.get(graphLabel).size(), graphLabel);
        }

        log.info("Starting server...");

        try {
            HttpServer server = HttpServer.create(new InetSocketAddress(PORT), 0);
            server.createContext("/",                new PageHandler("/html/index.html"));
            server.createContext("/d3.v3.min.js",    new PageHandler("/html/d3.v3.min.js"));
            server.createContext("/autocomplete.js", new PageHandler("/html/autocomplete.js"));
            server.createContext("/indiana.css",     new PageHandler("/html/indiana.css"));
            server.createContext("/contigs.csv",     new ContigsHandler());
            server.createContext("/contig",          new ContigHandler());
            server.createContext("/masked",          new MaskHandler());
            server.createContext("/edges",           new EdgeHandler());
            server.createContext("/links",           new LinksHandler());
            server.createContext("/link",            new LinkHandler());
            server.createContext("/dataset",         new DataSetHandler());
            server.setExecutor(null);
            server.start();
        } catch (IOException e) {
            throw new IndianaException("Unable to start server", e);
        }

        log.info("  listening on port {}", PORT);
    }
}
