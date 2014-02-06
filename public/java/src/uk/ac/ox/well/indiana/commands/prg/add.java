package uk.ac.ox.well.indiana.commands.prg;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;

import java.util.HashMap;

@Description(text="adds a FASTA file to the population reference graph")
public class add extends Module {
    @Argument(fullName="fasta", shortName="fa", doc="FASTA file")
    public HashMap<String, IndexedFastaSequenceFile> FASTAS;
    //public IndexedFastaSequenceFile FASTA;

    //@Argument(fullName="gff", shortName="g", doc="Gene feature format file")
    //public GFF3 GFF;

    @Override
    public void execute() {
        for (String key : FASTAS.keySet()) {
            log.info("{}={}", key, FASTAS.get(key));
        }
    }
}
