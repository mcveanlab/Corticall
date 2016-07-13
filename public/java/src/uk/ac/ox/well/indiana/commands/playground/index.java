package uk.ac.ox.well.indiana.commands.playground;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.File;
import java.io.PrintStream;

public class index extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM")
    public File SAMFILE;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
    }
}
