package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.util.ArrayList;
import java.util.List;

public class EvaluateSimContigs extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="links", shortName="l", doc="Links")
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(ROI.getColor(0).getSampleName());

        TraversalEngine el = new TraversalEngineFactory()
                .traversalColor(childColor)
                .graph(GRAPH)
                .links(LINKS)
                .make();

        TraversalEngine er = new TraversalEngineFactory()
                .traversalColor(childColor)
                .graph(GRAPH)
                .make();

        for (CortexRecord rr : ROI) {
            List<CortexVertex> elw = el.walk(rr.getKmerAsString());
            List<CortexVertex> erw = er.walk(rr.getKmerAsString());

            log.info("{} {} {}", rr.getKmerAsString(), elw.size(), erw.size());
        }
    }
}
