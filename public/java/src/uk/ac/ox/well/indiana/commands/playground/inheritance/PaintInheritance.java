package uk.ac.ox.well.indiana.commands.playground.inheritance;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;

import java.io.PrintStream;
import java.util.Arrays;

public class PaintInheritance extends Module {
    @Argument(fullName="reference", shortName="r", doc="Ref")
    public FastaSequenceFile REF;

    @Argument(fullName="mother", shortName="m", doc="Mother")
    public CortexGraph MOTHER;

    @Argument(fullName="father", shortName="f", doc="Father")
    public CortexGraph FATHER;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int kmerSize = MOTHER.getKmerSize();

        ReferenceSequence rs;
        while ((rs = REF.nextSequence()) != null) {
            String seq = rs.getBaseString();
            char[] painting = new char[seq.length() - kmerSize + 1];

            ProgressMeter pm = new ProgressMeter(log, 0, 10000, seq.length() - kmerSize, "Examining inheritance", "bases", " - ");
            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                pm.update();

                CortexKmer ck = new CortexKmer(seq.substring(i, i + kmerSize));

                boolean inMother = MOTHER.findRecord(ck) == null;
                boolean inFather = FATHER.findRecord(ck) == null;

                if      ( inMother && !inFather) { painting[i] = 'M'; }
                else if (!inMother &&  inFather) { painting[i] = 'F'; }
                else if (!inMother && !inFather) { painting[i] = 'N'; }
                else if ( inMother &&  inFather) { painting[i] = 'B'; }
            }

            out.println(rs.getName());
            out.println(new String(painting));
        }
    }
}
