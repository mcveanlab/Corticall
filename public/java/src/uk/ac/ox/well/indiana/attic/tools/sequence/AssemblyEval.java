package uk.ac.ox.well.indiana.attic.tools.sequence;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;

public class AssemblyEval extends Module {
    @Argument(fullName="fasta", shortName="f", doc="Fasta file to evaluate")
    public FastaSequenceFile FASTA;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Collection<String> seqs = new ArrayList<String>();

        ReferenceSequence seq;
        while ((seq = FASTA.nextSequence()) != null) {
            String sseq = new String(seq.getBases());

            seqs.add(sseq);
        }

        int numContigs = seqs.size();
        int n50Value = SequenceUtils.computeN50Value(seqs);

        out.println("numContigs: " + numContigs);
        out.println("n50: " + n50Value);
    }
}
