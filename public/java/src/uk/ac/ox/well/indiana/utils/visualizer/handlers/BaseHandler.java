package uk.ac.ox.well.indiana.utils.visualizer.handlers;

import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;

import java.io.IOException;
import java.io.OutputStream;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by kiran on 19/05/2017.
 */
public abstract class BaseHandler implements HttpHandler {
    public Map<String, String> query(String query) {
        Map<String, String> result = new HashMap<>();

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

