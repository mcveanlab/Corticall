package uk.ac.ox.well.indiana.attic.simulate;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.StringUtil;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
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

    @Argument(fullName="readLength", shortName="rl", doc="Read length to require")
    public Integer READ_LENGTH = 76;

    @Output(fullName="out_end1", shortName="o1", doc="File to which end1 should be written")
    public PrintStream o1;

    @Output(fullName="out_end2", shortName="o2", doc="File to which end2 should be written")
    public PrintStream o2;

    @Output
    public PrintStream out;

    private class PairedEndRead {
        private String end1;
        private String end2;
    }

    @Override
    public void execute() {
        Map<String, PairedEndRead> peReads = new HashMap<String, PairedEndRead>();

        int allReads = 0;
        int seenReads = 0;
        int writtenReads = 0;

        Map<Integer, String> qualStrings = new HashMap<Integer, String>();

        log.info("Processing reads...");
        for (SAMRecord read : BAM) {
            if (!read.getNotPrimaryAlignmentFlag() && read.getReadLength() == READ_LENGTH && read.getMappingQuality() > 0) {
                try {
                    String seq = new String(REF.getSubsequenceAt(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd()).getBases());

                    if (seq.length() == READ_LENGTH) {
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
                            o1.println(qualStrings.get(peReads.get(seqName).end1.length()));

                            o2.println("@" + seqName);
                            o2.println(peReads.get(seqName).end2);
                            o2.println("+");
                            o2.println(qualStrings.get(peReads.get(seqName).end2.length()));

                            peReads.remove(seqName);

                            writtenReads += 2;

                            if (N > 0 && writtenReads >= N) {
                                break;
                            }
                        }
                    }
                } catch (Exception e) {
                }

                seenReads++;
            }

            if (allReads % 1000000 == 0) {
                log.info("  processed {} reads", allReads);
            }
            allReads++;
        }

        TableWriter tw = new TableWriter(out);

        Map<String, String> te = new HashMap<String, String>();
        te.put("processed", String.valueOf(allReads));
        te.put("examined", String.valueOf(seenReads));
        te.put("wrote", String.valueOf(writtenReads));

        tw.addEntry(te);
    }
}
