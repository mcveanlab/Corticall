package uk.ac.ox.well.indiana.tools.examples;

import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;

import java.io.File;

public class HelloWorld extends Tool {
    @Argument(shortName="t", fullName="test", doc="A test argument")
    public String TEST = "value";

    @Argument(shortName="ra", fullName="requiredOpt", doc="Required option")
    public File REQUIRED_OPTION = null;

    @Argument(shortName="b", fullName="boolean", doc="flag")
    public Boolean FLAG = false;

    public int execute() {
        System.out.println(TEST);
        System.out.println(REQUIRED_OPTION);
        System.out.println(FLAG);

        return 0;
    }
}
