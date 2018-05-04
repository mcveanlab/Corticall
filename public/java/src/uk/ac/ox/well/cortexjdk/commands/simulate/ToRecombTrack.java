package uk.ac.ox.well.cortexjdk.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;

import java.io.PrintStream;

public class ToRecombTrack extends Module {
    @Argument(fullName="fasta", shortName="r", doc="Fasta")
    public FastaSequenceFile FASTA;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        ReferenceSequence rseq;
        while ((rseq = FASTA.nextSequence()) != null) {
            String chrName = rseq.getName().split(" ")[0];
            String seq = rseq.getBaseString();

            int start = 0;
            boolean lastHaplo = Character.isLowerCase(seq.charAt(start));
            for (int stop = 1; stop < seq.length(); stop++) {
                boolean thisHaplo = Character.isLowerCase(seq.charAt(stop));

                if (lastHaplo != thisHaplo || stop == seq.length() - 1) {
                    out.println(Joiner.on("\t").join(chrName, start, stop - 1, lastHaplo ? 0 : 1));

                    lastHaplo = thisHaplo;
                    start = stop;
                }
            }
        }
    }
}
