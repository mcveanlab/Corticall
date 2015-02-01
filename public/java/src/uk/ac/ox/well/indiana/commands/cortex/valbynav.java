package uk.ac.ox.well.indiana.commands.cortex;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;

import java.io.PrintStream;
import java.util.ArrayList;

@Description(text="Validates contigs by navigating the available Cortex graphs")
public class valbynav extends Module {
    @Argument(fullName="graph", shortName="g", doc="Cortex graph")
    public ArrayList<CortexGraph> GRAPHS;

    @Argument(fullName="contigs", shortName="c", doc="Contigs")
    public FastaSequenceFile CONTIGS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int kmerSize = GRAPHS.iterator().next().getKmerSize();

        ReferenceSequence rseq;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            String seq = new String(rseq.getBases());

            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                String sk = seq.substring(i, i + kmerSize);
                CortexKmer ck = new CortexKmer(sk);


            }
        }
    }
}
