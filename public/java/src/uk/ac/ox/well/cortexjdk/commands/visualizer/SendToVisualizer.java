package uk.ac.ox.well.cortexjdk.commands.visualizer;

import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.json.JSONObject;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.VisualizationStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;
import uk.ac.ox.well.cortexjdk.utils.visualizer.GraphVisualizer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

/**
 * Created by kiran on 15/09/2017.
 */
public class SendToVisualizer extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="kmer", shortName="k", doc="Kmer")
    public String KMER;

    @Argument(fullName="traversalColor", shortName="t", doc="Traversal color")
    public Integer TRAVERSAL_COLOR;

    @Argument(fullName="name", shortName="n", doc="Name to give to subgraph")
    public String NAME;

    @Argument(fullName="port", shortName="p", doc="Port")
    public Integer PORT = 9000;

    @Output
    public File out;

    @Override
    public void execute() {
        Collection<Integer> secondaryColors = new ArrayList<>();
        for (int c = 0; c < GRAPH.getNumColors(); c++) {
            if (c != TRAVERSAL_COLOR) {
                secondaryColors.add(c);
            }
        }

        TraversalEngine e = new TraversalEngineFactory()
                .graph(GRAPH)
                .traversalColor(TRAVERSAL_COLOR)
                .secondaryColors(secondaryColors)
                .traversalDirection(BOTH)
                .combinationOperator(OR)
                .stoppingRule(VisualizationStopper.class)
                .make();

        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = e.dfs(KMER);

        GraphVisualizer gv = new GraphVisualizer(PORT);
        JSONObject jo = gv.display(g, NAME);

        try {
            PrintStream jout = new PrintStream(out);

            jout.println(jo.toString(4));
            jout.flush();
        } catch (FileNotFoundException e1) {
            throw new CortexJDKException("File not found", e1);
        }
    }
}
