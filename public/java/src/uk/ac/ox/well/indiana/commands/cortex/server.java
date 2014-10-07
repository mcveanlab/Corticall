package uk.ac.ox.well.indiana.commands.cortex;

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
            log.info("Static request: {}", httpExchange);

            InputStream is = this.getClass().getResourceAsStream(page);

            BufferedReader in = new BufferedReader(new InputStreamReader(is));
            StringBuilder responseData = new StringBuilder();

            String line = null;
            while((line = in.readLine()) != null) {
                responseData.append(line).append("\n");
            }

            String response = responseData.toString();
            response = response.replaceAll("\\$NUM_CONTIGS", String.valueOf(contigs.size()));

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();
        }
    }

    private class JsonRequestHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {

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
            server.createContext("/test", new JsonRequestHandler());
            server.setExecutor(null);
            server.start();
        } catch (IOException e) {
            throw new IndianaException("Unable to start server", e);
        }

        log.info("  listening on port {}", PORT);
    }
}
