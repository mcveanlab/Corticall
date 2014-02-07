package uk.ac.ox.well.indiana.commands.prg;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

@Description(text="Extracts regions of interest from a FASTA file")
public class roi extends Module {
    @Argument(fullName="fasta", shortName="f", doc="FASTA file")
    public IndexedFastaSequenceFile FASTA;

    @Argument(fullName="gff", shortName="g", doc="GFF file")
    public GFF3 GFF;

    @Argument(fullName="ids", shortName="i", doc="Gene IDs to extract from GFF file", required=false)
    public HashSet<String> IDS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        for (GFF3Record gr : GFF) {
            if ("gene".equalsIgnoreCase(gr.getType())) {
                String id = gr.getAttribute("ID");

                if (IDS == null || IDS.contains(id)) {
                    String gene = SequenceUtils.extractGeneSequence(gr, FASTA);

                    out.println(">" + id + ".gene");
                    out.println(gene);

                    Collection<GFF3Record> exons = GFF3.getType("exon", GFF.getChildren(gr));
                    if (!exons.isEmpty()) {
                        String cds = SequenceUtils.extractCodingSequence(exons, FASTA);
                        String tr  = SequenceUtils.translateCodingSequence(cds);

                        out.println(">" + id + ".cds");
                        out.println(cds);

                        out.println(">" + id + ".tr");
                        out.println(tr);
                    }
                }
            }
        }
    }
}
