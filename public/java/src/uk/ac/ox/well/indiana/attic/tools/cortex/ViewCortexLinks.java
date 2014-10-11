package uk.ac.ox.well.indiana.attic.tools.cortex;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksRecord;

import java.io.PrintStream;

public class ViewCortexLinks extends Module {
    @Argument(fullName="cortexLinks", shortName="cl", doc="Cortex links")
    public CortexLinks CORTEX_LINKS;

    @Argument(fullName="headerOnly", shortName="H", doc="Only print the header")
    public Boolean HEADER_ONLY = false;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        if (HEADER_ONLY) {
            out.println(CORTEX_LINKS);
        } else {
            for (CortexLinksRecord cpr : CORTEX_LINKS) {
                out.println(cpr);
            }
        }
    }
}
