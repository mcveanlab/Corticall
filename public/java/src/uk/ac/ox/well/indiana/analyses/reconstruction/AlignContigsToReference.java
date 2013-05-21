package uk.ac.ox.well.indiana.analyses.reconstruction;

import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;

import java.io.File;

public class AlignContigsToReference extends Tool {
    @Argument(fullName="contigTable", shortName="ct", doc="Contig table")
    public File CONTIG_TABLE;

    //@Argument(fullName="")

    @Override
    public void execute() {
    }
}
