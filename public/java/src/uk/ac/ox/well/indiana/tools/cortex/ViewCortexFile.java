package uk.ac.ox.well.indiana.tools.cortex;

import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;

import java.io.PrintStream;

public class ViewCortexFile extends Tool {
    @Argument(fullName="cortexFile", shortName="c", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="headerOnly", shortName="H", doc="Only print the header")
    public Boolean HEADER_ONLY = false;

    @Output
    public PrintStream out;

    @Override
    public int execute() {
        if (HEADER_ONLY) {
            out.println(CORTEX_GRAPH);
        } else {
            for (CortexRecord cr : CORTEX_GRAPH) {
                out.println(cr);
            }
        }

        return 0;
    }
}
