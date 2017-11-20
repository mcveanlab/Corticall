package uk.ac.ox.well.cortexjdk.commands.utils;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;

import java.io.PrintStream;

public class Head extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="numRecords", shortName="n", doc="Num records")
    public Long NUM_RECORDS = 10L;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        for (long i = 0; i < (NUM_RECORDS < GRAPH.getNumRecords() ? NUM_RECORDS : GRAPH.getNumRecords()); i++) {
            out.println(GRAPH.getRecord(i));
        }
    }
}
