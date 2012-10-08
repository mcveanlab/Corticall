package uk.ac.ox.well.indiana.tools.cortex;

import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;

public class ViewCortexFile extends Tool {
    @Argument(fullName="cg", shortName="c", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Override
    public int execute() {
        for (CortexRecord cr : CORTEX_GRAPH) {
            System.out.println(cr);
        }

        return 0;
    }
}
