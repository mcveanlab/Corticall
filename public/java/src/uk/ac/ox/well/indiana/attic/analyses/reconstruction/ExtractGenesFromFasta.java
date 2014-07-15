package uk.ac.ox.well.indiana.attic.analyses.reconstruction;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;

import java.io.PrintStream;

public class ExtractGenesFromFasta extends Module {
    @Argument(fullName="reference", shortName="R", doc="Reference sequence")
    public IndexedFastaSequenceFile REFERENCE;

    @Argument(fullName="gff", shortName="gff", doc="GFF file")
    public GFF3 GFF;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        for (GFF3Record gr : GFF) {
            if ("gene".equals(gr.getType())) {
                String seq = new String(REFERENCE.getSubsequenceAt(gr.getSeqid(), gr.getStart(), gr.getEnd()).getBases());

                out.println(">" + gr.getAttribute("ID"));
                out.println(seq);
            }
        }
    }
}
