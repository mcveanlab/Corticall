package uk.ac.ox.well.indiana.attic.simulate;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

public class ExtractPerfectlyAligningReads extends Module {
    @Argument(fullName="sam", shortName="s", doc="SAM file")
    public SAMFileReader SAM;

    @Output(fullName="out_end1", shortName="o1", doc="File to which end1 should be written")
    public PrintStream o1;

    @Output(fullName="out_end2", shortName="o2", doc="File to which end2 should be written")
    public PrintStream o2;

    @Output
    public PrintStream out;

    private class AlignedPairedEndRead {
        public SAMRecord end1;
        public SAMRecord end2;
    }

    private boolean isPerfectRead(SAMRecord read) {
        String md = read.getStringAttribute("MD");

        return (md.equals(String.valueOf(read.getReadLength())));
    }

    @Override
    public void execute() {
        Map<String, AlignedPairedEndRead> peReads = new HashMap<String, AlignedPairedEndRead>();

        int allReads = 0;
        int perfectReads = 0;

        for (SAMRecord read : SAM) {
            String readName = read.getReadName();

            if (!peReads.containsKey(readName)) {
                peReads.put(readName, new AlignedPairedEndRead());

                peReads.get(readName).end1 = read;
            } else {
                peReads.get(readName).end2 = read;

                SAMRecord read1 = peReads.get(readName).end1;
                SAMRecord read2 = peReads.get(readName).end2;

                if (isPerfectRead(read1) && isPerfectRead(read2)) {
                    String seq1 = read1.getReadString();
                    String seq2 = read2.getReadString();

                    if (read1.getReadNegativeStrandFlag()) {
                        seq1 = SequenceUtils.reverseComplement(seq1);
                    }

                    if (read2.getReadNegativeStrandFlag()) {
                        seq2 = SequenceUtils.reverseComplement(seq2);
                    }

                    o1.println("@" + readName);
                    o1.println(seq1);
                    o1.println("+");
                    o1.println(read1.getBaseQualityString());

                    o2.println("@" + readName);
                    o2.println(seq2);
                    o2.println("+");
                    o2.println(read2.getBaseQualityString());

                    perfectReads++;
                }

                peReads.remove(readName);
            }

            allReads++;
        }

        TableWriter tw = new TableWriter(out);

        Map<String, String> te = new HashMap<String, String>();
        te.put("allReads", String.valueOf(allReads));
        te.put("perfectReads", String.valueOf(perfectReads));

        tw.addEntry(te);
    }
}
