package uk.ac.ox.well.indiana.commands.evaluate;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;

import java.io.PrintStream;

public class ExtractExons extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference")
    public IndexedFastaSequenceFile REF;

    @Argument(fullName="gff", shortName="g", doc="GFF")
    public GFF3 GFF;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        for (GFF3Record gr : GFF) {
            if (gr.getType().equals("exon")) {
                String id = gr.getAttribute("ID");
                String exon = new String(REF.getSubsequenceAt(gr.getSeqid(), gr.getStart(), gr.getEnd()).getBases());

                out.println(">" + id);
                out.println(exon);
            }
        }
    }
}
