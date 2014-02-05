package uk.ac.ox.well.indiana.attic.mod.qc;

import net.sf.samtools.SAMFileReader;
import uk.ac.ox.well.indiana.commands.Module;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: kiran
 * Date: 30/10/2013
 * Time: 18:36
 * To change this template use File | Settings | File Templates.
 */
public class NumRecords extends Module {
    @Override
    public void execute() {
        SAMFileReader sfr = new SAMFileReader(new File("test.bam"));
        //sfr.getBrowseableIndex().

        //To change body of implemented methods use File | Settings | File Templates.
    }
}
