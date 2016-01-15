package uk.ac.ox.well.indiana.commands.pacbio;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.PrintStream;

public class SliceReference extends Module {
    @Argument(fullName="ref", shortName="r", doc="Reference")
    public FastaSequenceFile REF;

    @Argument(fullName="length", shortName="l", doc="Length")
    public Integer LENGTH = 100;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        ReferenceSequence rseq;
        while ((rseq = REF.nextSequence()) != null) {
            String seq = new String(rseq.getBases());

            for (int i = 0; i + LENGTH < rseq.length(); i += LENGTH) {
                String slice = seq.substring(i, i + LENGTH);

                out.println(">" + rseq.getName() + "." + i);
                out.println(slice);
            }
        }
    }
}
