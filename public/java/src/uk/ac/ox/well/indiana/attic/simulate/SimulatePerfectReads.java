package uk.ac.ox.well.indiana.attic.simulate;

import net.sf.picard.PicardException;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

public class SimulatePerfectReads extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public SAMFileReader BAM;

    @Argument(fullName="reference", shortName="r", doc="Reference sequence")
    public IndexedFastaSequenceFile REF;

    @Output(fullName="out_end1", shortName="o1", doc="File to which end1 should be written")
    public PrintStream o1;

    @Output(fullName="out_end2", shortName="o2", doc="File to which end2 should be written")
    public PrintStream o2;

    private class PairedEndRead {
        private String end1;
        private String end2;
    }

    @Override
    public void execute() {
        Map<String, PairedEndRead> peReads = new HashMap<String, PairedEndRead>();

        int seenReads = 0;
        int writtenReads = 0;

        Map<Integer, String> qualStrings = new HashMap<Integer, String>();

        log.info("Processing reads...");
        for (SAMRecord read : BAM) {
            if (!qualStrings.containsKey(read.getReadLength())) {
                qualStrings.put(read.getReadLength(), StringUtil.repeatCharNTimes('I', read.getReadLength()));
            }

            String seqName = read.getReadName();

            try {
                String seq = new String(REF.getSubsequenceAt(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentStart() + read.getReadLength() - 1).getBases());

                if (!peReads.containsKey(seqName)) {
                    peReads.put(seqName, new PairedEndRead());

                    peReads.get(seqName).end1 = seq;
                } else {
                    peReads.get(seqName).end2 = seq;

                    o1.println("@" + seqName);
                    o1.println(peReads.get(seqName).end1);
                    o1.println("+");
                    o1.println(qualStrings.get(read.getReadLength()));

                    o2.println("@" + seqName);
                    o2.println(peReads.get(seqName).end2);
                    o2.println("+");
                    o2.println(qualStrings.get(read.getReadLength()));

                    peReads.remove(seqName);
                }

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
