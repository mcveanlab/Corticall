package uk.ac.ox.well.cortexjdk.utils.visualizer.handlers;

import com.google.api.client.http.HttpStatusCodes;
import com.sun.net.httpserver.Headers;
import com.sun.net.httpserver.HttpExchange;
import org.apache.commons.math3.util.Pair;

import java.io.*;
import java.nio.charset.Charset;

/**
 * Created by kiran on 19/05/2017.
 */
public class PageHandler extends BaseHandler {
    private boolean serveFromJar;

    public PageHandler(boolean serveFromFilesystem) {
        this.serveFromJar = serveFromFilesystem;
    }

    public Pair<Integer, String> respond(HttpExchange httpExchange) throws IOException {
        String page = "/html/" + httpExchange.getRequestURI();
        File f = new File("./public" + page);

        int code = HttpStatusCodes.STATUS_CODE_OK;
        String response;

        if (f.exists()) {
            InputStream is;

            if (!serveFromJar) {
                is = new FileInputStream(f);
            } else {
                is = this.getClass().getResourceAsStream(page);
            }

            BufferedReader in = new BufferedReader(new InputStreamReader(is, Charset.forName("UTF8")));
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
            code = HttpStatusCodes.STATUS_CODE_NOT_FOUND;
            response = "Page not found.";
        }

        return new Pair<>(code, response);
    }
}
