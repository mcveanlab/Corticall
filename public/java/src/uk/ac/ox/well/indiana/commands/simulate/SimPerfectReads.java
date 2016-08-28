package uk.ac.ox.well.indiana.commands.simulate;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.PrintStream;

public class SimPerfectReads extends Module {
    @Argument(fullName="ref", shortName="r", doc="Ref")
    public FastaSequenceFile REF;

    @Argument(fullName="readLength", shortName="l", doc="Read length")
    public Integer READ_LENGTH = 100;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        ReferenceSequence rs;
        while ((rs = REF.nextSequence()) != null) {
            String seq = new String(rs.getBases());

            for (int i = 0; i <= seq.length() - READ_LENGTH; i++) {
                out.println(">" + i);
                out.println(seq.substring(i, i + READ_LENGTH));
            }
        }
    }
}
