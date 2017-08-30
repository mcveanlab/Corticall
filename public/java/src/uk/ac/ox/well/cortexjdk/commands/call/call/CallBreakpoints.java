package uk.ac.ox.well.cortexjdk.commands.call.call;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.BubbleOpeningStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.AND;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

/**
 * Created by kiran on 30/08/2017.
 */
public class CallBreakpoints extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="links", shortName="l", doc="Links")
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="references", shortName="R", doc="References")
    public HashMap<String, KmerLookup> REFERENCES;

    @Output
    public PrintStream out;

    @Output(fullName="rout", shortName="ro", doc="Remaining ROI out")
    public File rout;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(new ArrayList<>(REFERENCES.keySet()));

        TraversalEngine eBreakpoint = new TraversalEngineFactory()
                .traversalColor(childColor)
                .joiningColors(parentColors)
                .traversalDirection(BOTH)
                .combinationOperator(AND)
                .stoppingRule(BubbleOpeningStopper.class)
                .graph(GRAPH)
                .links(LINKS)
                .rois(ROI)
                .make();
    }
}
