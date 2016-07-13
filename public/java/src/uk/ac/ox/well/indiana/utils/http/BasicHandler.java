package uk.ac.ox.well.indiana.utils.http;

import ch.qos.logback.classic.Logger;
import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.HashMap;
import java.util.Map;

public abstract class BasicHandler implements HttpHandler {
    private File documentRoot;
    private Logger log;

    public BasicHandler() {
        this.documentRoot = new File("./public/html/");
    }

    public BasicHandler(Logger log) {
        this.documentRoot = new File("./public/html/");
        this.log = log;
    }

    public BasicHandler(File documentRoot) {
        this.documentRoot = documentRoot;
        this.log = null;
    }

    public BasicHandler(File documentRoot, Logger log) {
        this.documentRoot = documentRoot;
        this.log = log;
    }

    private File getPage(HttpExchange httpExchange) {
        String page = httpExchange.getRequestURI().toString();

        if (page == null || page.isEmpty() || page.equals("/")) {
            page = "index.html";
        }

        return new File(documentRoot.getAbsolutePath() + "/" + page);
    }

    private Map<String, String> getQuery(HttpExchange httpExchange) {
        String query = httpExchange.getRequestURI().getQuery();

        Map<String, String> result = new HashMap<>();

        if (query != null) {
            for (String param : query.split("&")) {
                String pair[] = param.split("=", 2);
                if (pair.length > 1) {
                    result.put(pair[0], pair[1]);
                } else {
                    result.put(pair[0], "");
                }
            }
        }

        return result;
    }

    private void write(HttpExchange httpExchange, String response) throws IOException {
        int code = 200;
        if (response == null || response.isEmpty()) {
            code = 404;
            response = "Not found.";
        }

        httpExchange.sendResponseHeaders(code, response.length());

        OutputStream os = httpExchange.getResponseBody();
        os.write(response.getBytes());
        os.close();
    }

    @Override
    public void handle(HttpExchange httpExchange) throws IOException {
        File page = getPage(httpExchange);
        Map<String, String> query = getQuery(httpExchange);

        if (log != null) log.info("Request: page={} query={}", page, query);

        String result = process(page, query);

        write(httpExchange, result);
    }

    public abstract String process(File page, Map<String, String> query);
}

