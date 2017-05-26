package uk.ac.ox.well.indiana.commands.visualizer;

/**
 * Created by kiran on 10/05/2017.
 */
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.visualizer.GraphVisualizationFactory;
import uk.ac.ox.well.indiana.utils.visualizer.GraphVisualizer;

@Description(text="Starts the server for visualizing assembly data")
public class VisualCortex extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="rois", shortName="r", doc="ROIs")
    public CortexGraph ROIS;

    @Argument(fullName="port", shortName="p", doc="Port")
    public Integer PORT = 9000;

    @Override
    public void execute() {
        log.info("Graph: {} {}", GRAPH.getNumRecords(), GRAPH.getCortexFile());
        log.info("Rois:  {} {}", ROIS.getNumRecords(), ROIS.getCortexFile());
        log.info("");

        GraphVisualizer gv = new GraphVisualizationFactory()
                .port(PORT)
                .logger(log)
                .graph(GRAPH)
                .rois(ROIS)
                .make();

    }
}
