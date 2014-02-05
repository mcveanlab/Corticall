package uk.ac.ox.well.indiana.attic.analyses.roi;

import net.sf.picard.reference.FastaSequenceFile;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;

public class FindRelevantContigs extends Module {
    @Argument(fullName="sequences", shortName="s", doc="Sequences FASTA")
    public FastaSequenceFile FASTA;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Override
    public void execute() {
        //To change body of implemented methods use File | Settings | File Templates.
    }
}
