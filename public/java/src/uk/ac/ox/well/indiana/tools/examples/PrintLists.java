package uk.ac.ox.well.indiana.tools.examples;

import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;

public class PrintLists extends Tool {
    @Argument(fullName="list", shortName="l", doc="A list to load and print")
    public ArrayList<String> LIST;

    @Output
    public PrintStream out;

    @Override
    public int execute() {
        log.info("This is a PrintLists test");

//        for (CortexGraph entry : LIST) {
//            out.println(entry.getCortexFile().getName());
//        }

        for (String entry : LIST) {
            out.println(entry);
        }

        return 0;
    }
}
