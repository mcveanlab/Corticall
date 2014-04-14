package uk.ac.ox.well.indiana.attic.tools.sequence;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.PrintStream;
import java.util.ArrayList;

@Description(text="Removes contigs with empty sequences from FASTA file(s)")
public class RemoveEmptyFastaSequences extends Module {
    @Argument(fullName="fasta", shortName="f", doc="FASTA file")
    public ArrayList<FastaSequenceFile> FASTAS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        log.info("Processing FASTA files...");

        int index = 0;
        for (FastaSequenceFile f : FASTAS) {
            log.info("  FASTA {}/{}", index, FASTAS.size());
            index++;

            ReferenceSequence rseq;

            while ((rseq = f.nextSequence()) != null) {
                if (rseq.length() > 0) {
                    out.println(">" + rseq.getName());
                    out.println(new String(rseq.getBases()));
                }
            }
        }
    }
}
