package uk.ac.ox.well.indiana.attic.simulate;

import net.sf.picard.PicardException;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

public class SimulatePerfectReads extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public SAMFileReader BAM;

    @Argument(fullName="reference", shortName="r", doc="Reference sequence")
    public IndexedFastaSequenceFile REF;

    @Argument(fullName="numReads", shortName="n", doc="Number of reads to process")
    public Integer N = 10;

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
            if (!read.getNotPrimaryAlignmentFlag()) {
                try {
                    String seq = new String(REF.getSubsequenceAt(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd()).getBases());

                    if (seq.length() == read.getReadString().length()) {
                        if (read.getReadNegativeStrandFlag()) {
                            seq = SequenceUtils.reverseComplement(seq);
                        }

                        if (!qualStrings.containsKey(seq.length())) {
                            qualStrings.put(seq.length(), StringUtil.repeatCharNTimes('I', seq.length()));
                        }

                        String seqName = read.getReadName();
                        if (!peReads.containsKey(seqName)) {
                            peReads.put(seqName, new PairedEndRead());

                            peReads.get(seqName).end1 = seq;
                        } else {
                            peReads.get(seqName).end2 = seq;

                            o1.println("@" + seqName);
                            o1.println(peReads.get(seqName).end1);
                            o1.println("+");
                            o1.println(qualStrings.get(seq.length()));

                            o2.println("@" + seqName);
                            o2.println(peReads.get(seqName).end2);
                            o2.println("+");
                            o2.println(qualStrings.get(seq.length()));

                            peReads.remove(seqName);

                            writtenReads++;

                            if (N > 0 && writtenReads >= N) {
                                break;
                            }
                        }
                    }
                } catch (PicardException e) {
                    //log.warn("Discarded exception: {}", e);
                }
            }

            if (seenReads % 1000000 == 0) {
                log.info("  processed {} reads, wrote {} perfect reads", seenReads, writtenReads);
            }
            seenReads++;
        }

        log.info("Saw {} reads, wrote {} perfect reads.", seenReads, writtenReads);
    }
}
