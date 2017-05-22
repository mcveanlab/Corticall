package uk.ac.ox.well.indiana.utils.visualizer;

import org.jgrapht.DirectedGraph;
import org.slf4j.Logger;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;

/**
 * Created by kiran on 19/05/2017.
 */
public class VisualCortexConfiguration {
    private int port;
    private Logger log;
    private CortexGraph graph;
    private DirectedGraph<CortexVertex, CortexEdge> subGraph;

    public int getPort() { return port; }
    public void setPort(int port) { this.port = port; }

    public Logger getLogger() { return log; }
    public void setLogger(Logger log) { this.log = log; }

    public CortexGraph getGraph() { return graph; }
    public void setGraph(CortexGraph graph) { this.graph = graph; }

    public DirectedGraph<CortexVertex, CortexEdge> getSubGraph() { return subGraph; }
    public void setSubGraph(DirectedGraph<CortexVertex, CortexEdge> subGraph) { this.subGraph = subGraph; }
}
