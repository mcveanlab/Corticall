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
            int i = 0;
            for (CortexRecord cr : CORTEX_GRAPH) {
                i++;

                if (i % (CORTEX_GRAPH.getNumRecords()/10) == 0) {
                    log.info("{}/{} records processed", i, CORTEX_GRAPH.getNumRecords());
                }

                if (i > 6860970 - 10 && i < 6860970 + 10) {
                    log.info("{}: {}", i, cr);
                }

                /*
                if (satisfiesConstraints(cr)) {
                    String homeName = getKmerHomeContigName(cr);

                    if (homeName == null) {
                        out.println(cr);
                    } else {
                        out.println(homeName + " " + cr);
                    }
                }
                */
            }

            log.info("Loaded {}/{} records total", i, CORTEX_GRAPH.getNumRecords());
        }

        return 0;
    }
}
