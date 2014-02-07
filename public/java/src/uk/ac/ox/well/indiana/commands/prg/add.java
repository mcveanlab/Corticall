package uk.ac.ox.well.indiana.commands.prg;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;

import java.util.ArrayList;
import java.util.HashMap;

@Description(text="adds FASTA file(s) to the population reference graph")
public class add extends Module {
    @Argument(fullName="fasta", shortName="f", doc="FASTA file")
    public IndexedFastaSequenceFile FASTA;

    @Argument(fullName="gff", shortName="g", doc="GFF file")
    public GFF3 GFF;

    @Argument(fullName="genelist", shortName="l", doc="List of genes to include in graph")
    public ArrayList<String> GENES = new ArrayList<String>();

    @Override
    public void execute() {
    }
}
