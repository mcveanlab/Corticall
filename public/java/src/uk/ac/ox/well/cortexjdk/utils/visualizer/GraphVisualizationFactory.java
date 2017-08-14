package uk.ac.ox.well.cortexjdk.utils.visualizer;

import org.jgrapht.DirectedGraph;
import org.slf4j.Logger;
import uk.ac.ox.well.cortexjdk.Main;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;

/**
 * Created by kiran on 19/05/2017.
 */
public class GraphVisualizationFactory {
    private GraphVisualizationConfiguration vcc;

    public GraphVisualizationFactory() {
        vcc = new GraphVisualizationConfiguration();

        vcc.setPort(9000);
        vcc.setLogger(Main.getLogger());
    }

    public GraphVisualizationFactory port(int port) { vcc.setPort(port); return this; }
    public GraphVisualizationFactory logger(Logger log) { vcc.setLogger(log); return this; }
    public GraphVisualizationFactory graph(CortexGraph graph) { vcc.setGraph(graph); return this; }
    public GraphVisualizationFactory rois(CortexGraph rois) { vcc.setRois(rois); return this; }
    public GraphVisualizationFactory subgraph(DirectedGraph<CortexVertex, CortexEdge> g) { vcc.setSubGraph(g); return this; }

    public GraphVisualizer make() { return new GraphVisualizer(vcc); }
}
