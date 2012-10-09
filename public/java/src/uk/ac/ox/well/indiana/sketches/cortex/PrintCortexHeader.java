package uk.ac.ox.well.indiana.sketches.cortex;

import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;

import java.io.PrintStream;

public class PrintCortexHeader extends Tool {
    @Argument(fullName="cortexFile", shortName="cf", doc="Cortex file")
    public CortexGraph CORTEX_FILE;

    @Output
    public PrintStream out;

    @Override
    public int execute() {
        out.println(CORTEX_FILE);

        return 0;
    }
}
