package uk.ac.ox.well.indiana.commands.index;

import htsjdk.samtools.SAMFileReader;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;

public class ReadNameIndex extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public SAMFileReader BAM;

    @Override
    public void execute() {

    }
}
