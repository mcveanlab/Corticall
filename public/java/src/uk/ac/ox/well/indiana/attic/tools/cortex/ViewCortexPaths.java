package uk.ac.ox.well.indiana.attic.tools.cortex;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.paths.CortexPaths;
import uk.ac.ox.well.indiana.utils.io.cortex.paths.CortexPathsRecord;

import java.io.PrintStream;

public class ViewCortexPaths extends Module {
    @Argument(fullName="cortexPaths", shortName="cp", doc="Cortex paths")
    public CortexPaths CORTEX_PATHS;

    @Argument(fullName="headerOnly", shortName="H", doc="Only print the header")
    public Boolean HEADER_ONLY = false;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        if (HEADER_ONLY) {
            out.println(CORTEX_PATHS);
        } else {
            for (CortexPathsRecord cpr : CORTEX_PATHS) {
                out.println(cpr);
            }
        }
    }
}
