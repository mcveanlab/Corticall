package uk.ac.ox.well.indiana.tools.examples;

import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

public class PrintLists extends Tool {
    @Argument(fullName="list", shortName="l", doc="A list to load and print")
    public ArrayList<String> LIST;

    @Output
    public PrintStream out;

    @Override
    public int execute() {
        for (String entry : LIST) {
            System.out.println(entry);
        }

        return 0;
    }
}
