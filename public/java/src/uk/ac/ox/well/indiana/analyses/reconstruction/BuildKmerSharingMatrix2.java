package uk.ac.ox.well.indiana.analyses.reconstruction;

import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.File;
import java.io.PrintStream;

public class BuildKmerSharingMatrix2 extends Tool {
    @Argument(fullName="contigsTable", shortName="ct", doc="Contigs table")
    public File CONTIG_TABLE;

    @Argument(fullName="ignoreNonPanelKmers", shortName="inpk", doc="Only examine sharing for panel kmers.  Ignore non-panel kmers.")
    public Boolean IGNORE_NON_PANEL_KMERS = false;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
    }
}
