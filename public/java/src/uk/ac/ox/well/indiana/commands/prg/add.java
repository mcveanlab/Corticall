package uk.ac.ox.well.indiana.commands.prg;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;

import java.util.HashMap;

@Description(text="adds a FASTA file to the population reference graph")
public class add extends Module {
    @Argument(fullName="fasta", shortName="fa", doc="name:FASTA key-value pair")
    public HashMap<String, IndexedFastaSequenceFile> FASTAS;

    @Argument(fullName="gff", shortName="g", doc="name:GFF key-value pair", required=false)
    public HashMap<String, GFF3> GFFS;

    @Override
    public void execute() {
        for (String key : FASTAS.keySet()) {
            log.info("{}={}", key, FASTAS.get(key));
        }

        for (String key : GFFS.keySet()) {
            log.info("{}={}", key, GFFS.get(key));
        }
    }
}
