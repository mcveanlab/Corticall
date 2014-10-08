package uk.ac.ox.well.indiana.commands.cortex;

import com.sun.net.httpserver.Headers;
import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;

import java.io.*;
import java.net.InetSocketAddress;
import java.util.*;

@Description(text="Starts the server for visualizing assembly data")
public class server extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs (FASTA)")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="graph", shortName="g", doc="Cortex graph file")
    public HashMap<String, File> GRAPHS;

    @Argument(fullName="port", shortName="p", doc="Port")
    public Integer PORT = 9000;

    private Map<String, String> contigs;

    private class StaticPageHandler implements HttpHandler {
        private String page;

        public StaticPageHandler(String page) {
            this.page = page;
        }

        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            InputStream is = this.getClass().getResourceAsStream(page);

            BufferedReader in = new BufferedReader(new InputStreamReader(is));
            StringBuilder responseData = new StringBuilder();

            String line = null;
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

    private class ContigsRequestHandler implements HttpHandler {
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

    private class ContigRequestHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = queryToMap(httpExchange.getRequestURI().getQuery());

            StringBuilder rb = new StringBuilder();
            rb.append("{\"seq\": \"").append(contigs.get(query.get("contigName"))).append("\"}\n");

            String response = rb.toString();

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("GET contig seq    : {}, response length: {}", query.get("contigName"), response.length());
        }
    }

    private Map<String, String> loadContigs() {
        Map<String, String> contigs = new LinkedHashMap<String, String>();

        ReferenceSequence rseq;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            String name = rseq.getName();

            contigs.put(name, new String(rseq.getBases()));
        }

        return contigs;
    }

    @Override
    public void execute() {
        log.info("Loading contigs...");
        contigs = loadContigs();
        log.info("  loaded {} contigs", contigs.size());

        log.info("Starting server...");

        try {
            HttpServer server = HttpServer.create(new InetSocketAddress(PORT), 0);
            server.createContext("/", new StaticPageHandler("/html/index.html"));
            server.createContext("/d3.v3.min.js", new StaticPageHandler("/html/d3.v3.min.js"));
            server.createContext("/autocomplete.js", new StaticPageHandler("/html/autocomplete.js"));
            server.createContext("/contigs.csv", new ContigsRequestHandler());
            server.createContext("/contig", new ContigRequestHandler());
            server.setExecutor(null);
            server.start();
        } catch (IOException e) {
            throw new IndianaException("Unable to start server", e);
        }

        log.info("  listening on port {}", PORT);
    }
}
