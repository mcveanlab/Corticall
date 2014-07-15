package uk.ac.ox.well.indiana.attic.analyses.nahr;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.PrintStream;

public class TrimContigEnds extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs (FASTA)")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="trim", shortName="t", doc="Number of bases to trim from both the left and right sides of the contig")
    public Integer TRIM = 5;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        ReferenceSequence rseq;

        while ((rseq = CONTIGS.nextSequence()) != null) {
            String seq = new String(rseq.getBases());

            String subseq = seq.substring(TRIM, seq.length() - TRIM);

            out.println(">" + rseq.getName());
            out.println(subseq);
        }
    }
}
