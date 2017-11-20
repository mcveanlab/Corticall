package uk.ac.ox.well.cortexjdk.commands.visualizer;

/**
 * Created by kiran on 10/05/2017.
 */
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Description;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.visualizer.GraphVisualizationFactory;
import uk.ac.ox.well.cortexjdk.utils.visualizer.GraphVisualizer;

@Description(text="Starts the server for visualizing assembly data")
public class VisualCortex extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph", required=false)
    public CortexGraph GRAPH;

    @Argument(fullName="rois", shortName="r", doc="ROIs", required=false)
    public CortexGraph ROIS;

    @Argument(fullName="port", shortName="p", doc="Port")
    public Integer PORT = 9000;

    @Override
    public void execute() {
        GraphVisualizationFactory gvf = new GraphVisualizationFactory()
            .port(PORT)
            .logger(log);

        if (GRAPH != null) {
            log.info("Graph: {} {}", GRAPH.getNumRecords(), GRAPH.getFile());

            gvf = gvf.graph(GRAPH);
        }
        if (ROIS != null) {
            log.info("Rois:  {} {}", ROIS.getNumRecords(), ROIS.getFile());

            gvf = gvf.rois(ROIS);
        }

        GraphVisualizer gv = gvf.make();
    }
}
