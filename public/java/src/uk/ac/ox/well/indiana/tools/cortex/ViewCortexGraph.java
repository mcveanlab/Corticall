package uk.ac.ox.well.indiana.tools.cortex;

import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;

public class ViewCortexGraph extends ViewCortexBase {
    @Argument(fullName="headerOnly", shortName="H", doc="Only print the header")
    public Boolean HEADER_ONLY = false;

    @Override
    public int execute() {
        if (HEADER_ONLY) {
            out.println(CORTEX_GRAPH);
        } else {
            for (CortexRecord cr : CORTEX_GRAPH) {
                if (satisfiesConstraints(cr)) {
                    out.println(cr);
                }
            }
        }

        return 0;
    }
}
