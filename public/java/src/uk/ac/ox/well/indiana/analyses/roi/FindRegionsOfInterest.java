package uk.ac.ox.well.indiana.analyses.roi;

import net.sf.samtools.SAMFileReader;
import uk.ac.ox.well.indiana.tools.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;

import java.util.ArrayList;

public class FindRegionsOfInterest extends Module {
    @Argument(fullName="b", shortName="bam", doc="BAM file(s) to process")
    public ArrayList<SAMFileReader> BAMS;

    @Override
    public void execute() {

    }
}
