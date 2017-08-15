package uk.ac.ox.well.cortexjdk.commands.utils;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;

import java.io.PrintStream;

public class Range extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="start", shortName="s", doc="Start record")
    public Long START = 0L;

    @Argument(fullName="end", shortName="e", doc="End record")
    public Long END = 0L;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        for (long i = START; i < END; i++) {
            out.println(GRAPH.getRecord(i));
        }
    }
}
