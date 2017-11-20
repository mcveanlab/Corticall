package uk.ac.ox.well.cortexjdk.utils.visualizer;

import org.jgrapht.DirectedGraph;
import org.slf4j.Logger;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;

/**
 * Created by kiran on 19/05/2017.
 */
public class GraphVisualizationConfiguration {
    private int port;
    private Logger log;
    private CortexGraph graph;
    private CortexGraph rois;
    private DirectedGraph<CortexVertex, CortexEdge> subGraph;

    public int getPort() { return port; }
    public void setPort(int port) { this.port = port; }

    public Logger getLogger() { return log; }
    public void setLogger(Logger log) { this.log = log; }

    public CortexGraph getGraph() { return graph; }
    public void setGraph(CortexGraph graph) { this.graph = graph; }

    public DirectedGraph<CortexVertex, CortexEdge> getSubGraph() { return subGraph; }
    public void setSubGraph(DirectedGraph<CortexVertex, CortexEdge> subGraph) { this.subGraph = subGraph; }

    public CortexGraph getRois() { return rois; }
    public void setRois(CortexGraph rois) { this.rois = rois; }
}
