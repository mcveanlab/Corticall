package uk.ac.ox.well.indiana.attic.analyses.nahr;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;

public class TraverseIncorrectly extends Module {
    @Argument(fullName="seq", shortName="s", doc="Sequence considered correct")
    public String SEQUENCE;

    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Override
    public void execute() {

    }
}
