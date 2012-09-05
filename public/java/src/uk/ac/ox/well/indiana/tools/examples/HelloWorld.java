package uk.ac.ox.well.indiana.tools.examples;

import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;

public class HelloWorld extends Tool {
    @Argument(shortName="t", fullName="test", doc="A test argument")
    public String TEST = "value";

    public int execute() {
        System.out.println(TEST);

        return 0;
    }
}
