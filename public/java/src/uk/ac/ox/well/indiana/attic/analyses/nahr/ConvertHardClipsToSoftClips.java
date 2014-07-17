package uk.ac.ox.well.indiana.attic.analyses.nahr;

import htsjdk.samtools.*;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class ConvertHardClipsToSoftClips extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public ArrayList<SAMFileReader> BAMS;

    @Output
    public File out;

    @Override
    public void execute() {
        Set<SAMFileHeader> headers = new HashSet<SAMFileHeader>();
        for (SAMFileReader bam : BAMS) {
            headers.add(bam.getFileHeader());
        }

        SamFileHeaderMerger sfhm = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate, headers, true);
        SAMFileWriterFactory sfwf = new SAMFileWriterFactory().setCreateIndex(true);
        SAMFileWriter sfw = sfwf.makeBAMWriter(sfhm.getMergedHeader(), false, out);

        log.info("Processing reads...");
        int changedRecords = 0;
        int totalRecords = 0;
        for (SAMFileReader bam : BAMS) {
            for (SAMRecord read : bam) {
                List<CigarElement> ces = new ArrayList<CigarElement>();

                boolean changed = false;

                for (CigarElement ce : read.getCigar().getCigarElements()) {
                    if (ce.getOperator().equals(CigarOperator.H)) {
                        ces.add(new CigarElement(ce.getLength(), CigarOperator.S));

                        changed = true;
                    } else {
                        ces.add(ce);
                    }
                }

                read.setCigar(new Cigar(ces));

                sfw.addAlignment(read);

                if (changed) { changedRecords++; }
                totalRecords++;
            }
        }

        log.info("  changed {}/{} records", changedRecords, totalRecords);

        sfw.close();
    }
}
