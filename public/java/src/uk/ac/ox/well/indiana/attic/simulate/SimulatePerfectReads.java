package uk.ac.ox.well.indiana.attic.simulate;

import net.sf.picard.PicardException;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.PrintStream;

public class SimulatePerfectReads extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public SAMFileReader BAM;

    @Argument(fullName="reference", shortName="r", doc="Reference sequence")
    public IndexedFastaSequenceFile REF;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int seenReads = 0;
        int writtenReads = 0;

        log.info("Processing reads...");
        for (SAMRecord read : BAM) {
            String seqName = read.getReadName();

            try {
                String seq = new String(REF.getSubsequenceAt(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentStart() + read.getReadLength() - 1).getBases());

                out.println(">" + seqName);
                out.println(seq);

                writtenReads++;
            } catch (PicardException e) {
                //log.warn("Discarded exception: {}", e);
            }

            if (seenReads % 1000000 == 0) {
                log.info("  processed {} reads", seenReads);
            }
            seenReads++;
        }

        log.info("Saw {} reads, wrote {} perfect reads.", seenReads, writtenReads);
    }
}
