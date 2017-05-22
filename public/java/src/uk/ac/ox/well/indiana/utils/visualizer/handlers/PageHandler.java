package uk.ac.ox.well.indiana.utils.visualizer.handlers;

import com.sun.net.httpserver.Headers;
import com.sun.net.httpserver.HttpExchange;

import java.io.*;

/**
 * Created by kiran on 19/05/2017.
 */
public class PageHandler extends BaseHandler {
    private boolean serveFromJar;

    public PageHandler(boolean serveFromFilesystem) {
        this.serveFromJar = serveFromFilesystem;
    }

    @Override
    public void handle(HttpExchange httpExchange) throws IOException {
        String page = "/html/" + httpExchange.getRequestURI();
        File f = new File("./public" + page);

        int code = 200;
        String response;

        if (f.exists()) {
            InputStream is;

            if (!serveFromJar) {
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
