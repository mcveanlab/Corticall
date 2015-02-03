package uk.ac.ox.well.indiana.commands.simulate;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.PrintStream;

public class CreateSequenceDictionary extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs")
    public FastaSequenceFile CONTIGS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        out.println("@HD\tVN:1.0\tSO:unsorted");

        ReferenceSequence rseq;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            out.println("@SQ\tSN:" + rseq.getName() + "\tLN:" + rseq.length());
        }
    }
}
