package uk.ac.ox.well.indiana.utils.visualizer;

import org.jgrapht.DirectedGraph;
import org.slf4j.Logger;
import uk.ac.ox.well.indiana.Main;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;

/**
 * Created by kiran on 19/05/2017.
 */
public class VisualCortexFactory {
    private VisualCortexConfiguration vcc;

    public VisualCortexFactory() {
        vcc = new VisualCortexConfiguration();

        vcc.setPort(9000);
        vcc.setLogger(Main.getLogger());
    }

    public VisualCortexFactory port(int port) { vcc.setPort(port); return this; }
    public VisualCortexFactory logger(Logger log) { vcc.setLogger(log); return this; }
    public VisualCortexFactory graph(CortexGraph graph) { vcc.setGraph(graph); return this; }
    public VisualCortexFactory subgraph(DirectedGraph<CortexVertex, CortexEdge> g) { vcc.setSubGraph(g); return this; }

    public VisualCortex make() { return new VisualCortex(vcc); }
}
