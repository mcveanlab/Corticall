package uk.ac.ox.well.indiana.tools.cortex;

import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;

import java.io.PrintStream;
import java.util.Date;
import java.util.TreeSet;

public class SortCortexRecords extends Tool {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Output
    public PrintStream out;

    @Override
    public int execute() {
        TreeSet<CortexRecord> sortedRecords = new TreeSet<CortexRecord>();

        long recordsSeen = 0;
        long recordsWritten = 0;
        long recordsTotal = CORTEX_GRAPH.getNumRecords();

        for (CortexRecord cr : CORTEX_GRAPH) {
            sortedRecords.add(cr);

            recordsSeen++;
            if (recordsSeen % 100000 == 0) {
                System.out.println("records seen: " + recordsSeen + "/" + recordsTotal);
            }
        }

        System.out.println("loaded " + recordsSeen + "/" + recordsTotal);

        out.println(CORTEX_GRAPH);

        for (CortexRecord cr : sortedRecords) {
            out.println(cr);

            recordsWritten++;
            if (recordsWritten % 100000 == 0) {
                System.out.println("records written: " + recordsWritten + "/" + recordsTotal);
            }
        }

        System.out.println("wrote " + recordsWritten + "/" + recordsTotal);

        return 0;
    }
}
