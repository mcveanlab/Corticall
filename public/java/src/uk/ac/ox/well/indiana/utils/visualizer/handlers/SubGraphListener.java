package uk.ac.ox.well.indiana.utils.visualizer.handlers;

import com.google.api.client.http.HttpStatusCodes;
import com.sun.net.httpserver.HttpExchange;
import org.apache.commons.math3.util.Pair;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * Created by kiran on 29/05/2017.
 */
public class SubGraphListener extends BaseHandler {
    private String preLoadedGraph = "";

    @Override
    public Pair<Integer, String> respond(HttpExchange httpExchange) throws IOException {
        if (httpExchange.getRequestMethod().equals("POST")) {
            StringBuilder sb = new StringBuilder();

            new BufferedReader(new InputStreamReader(httpExchange.getRequestBody()))
                    .lines()
                    .forEach((String s) -> sb.append(s + "\n"));

            preLoadedGraph = sb.toString();
        }

        return new Pair<>(HttpStatusCodes.STATUS_CODE_OK, preLoadedGraph);
    }
}
