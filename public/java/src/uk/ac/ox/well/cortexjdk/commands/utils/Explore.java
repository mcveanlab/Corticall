package uk.ac.ox.well.cortexjdk.commands.utils;

import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.DestinationStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.*;

import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

public class Explore extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="links", shortName="l", doc="Links")
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName="sample", shortName="s", doc="Sample")
    public TreeSet<String> SAMPLES;

    @Argument(fullName="begin", shortName="b", doc="Beginning kmer")
    public String BEGIN;

    @Argument(fullName="end", shortName="e", doc="Ending kmer")
    public String END;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TraversalEngine e = new TraversalEngineFactory()
                .traversalColors(GRAPH.getColorForSampleName("3D7"))
                .traversalDirection(BOTH)
                .combinationOperator(OR)
                .stoppingRule(DestinationStopper.class)
                .graph(GRAPH)
                .links(LINKS)
                .make();

        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = e.dfs(BEGIN, END);
        if (g == null) {
            g = e.dfs(END, BEGIN);
        }

        List<CortexVertex> w = TraversalUtils.toWalk(g, BEGIN, GRAPH.getColorForSampleName("3D7"));

        /*
        Map<Integer, List<List<List<CortexVertex>>>> m = TraversalUtils.fill(w, GRAPH, BACKGROUND_LINKS, new HashSet<>(GRAPH.getColorsForSampleNames(SAMPLES)), 0);

        for (int c : m.keySet()) {
            String sample = GRAPH.getSampleName(c);

            List<List<List<CortexVertex>>> lllc = m.get(c);

            for (int i = 0; i < lllc.size(); i++) {
                for (int j = 0; j < lllc.get(i).size(); j++) {
                    log.info("{} {} {} {} {}", String.format("%-2d", c), String.format("%-20s", sample), i, j, TraversalUtils.toContig(lllc.get(i).get(j)));
                }
            }

            log.info("--");
        }
        */
    }
}
