package uk.ac.ox.well.indiana.utils.visualizer;

import com.sun.net.httpserver.HttpServer;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.visualizer.handlers.PageHandler;
import uk.ac.ox.well.indiana.utils.visualizer.handlers.SubGraphHandler;

import java.io.*;
import java.net.InetSocketAddress;

public class GraphVisualizer {
    private GraphVisualizationConfiguration vcc;

    public GraphVisualizer(GraphVisualizationConfiguration vcc) {
        this.vcc = vcc;

        try {
            HttpServer server = HttpServer.create(new InetSocketAddress(vcc.getPort()), 0);
            server.createContext("/", new PageHandler(false));
            server.createContext("/graph", new SubGraphHandler(vcc.getGraph(), vcc.getRois()));
            server.setExecutor(null);
            server.start();
        } catch (IOException e) {
            throw new IndianaException("Unable to start server", e);
        }

        vcc.getLogger().info("Started server on port {}...", vcc.getPort());
    }
}
