package uk.ac.ox.well.indiana.attic.analyses.roi;

import htsjdk.samtools.SAMFileReader;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;

import java.util.ArrayList;

public class FindRegionsOfInterest extends Module {
    @Argument(fullName="b", shortName="bam", doc="BAM file(s) to process")
    public ArrayList<SAMFileReader> BAMS;

    @Override
    public void execute() {

    }
}
