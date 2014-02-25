package uk.ac.ox.well.indiana.attic.analyses.LongContigs;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.PrintStream;

public class KmerizeReference extends Module {
    @Argument(fullName="fasta", shortName="f", doc="Reference genome")
    public FastaSequenceFile FASTA;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        ReferenceSequence rseq;
        while ((rseq = FASTA.nextSequence()) != null) {
            String[] name = rseq.getName().split("\\s+");
            String seq = new String(rseq.getBases());

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                String kmer = seq.substring(i, i + KMER_SIZE);

                out.println(">" + name[0] + "." + i);
                out.println(kmer);
            }
        }
    }
}
