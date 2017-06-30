package uk.ac.ox.well.indiana.commands.playground.assemblies;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;

/**
 * Created by kiran on 30/06/2017.
 */
public class ExtractExons extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference")
    public IndexedFastaSequenceFile REFERENCE;

    @Argument(fullName="gff", shortName="g", doc="GFF3 file")
    public GFF3 GFF;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        for (GFF3Record gr : GFF) {
            if (gr.getType().equals("exon")) {
                String id = gr.getAttribute("ID");
                String seq = REFERENCE.getSubsequenceAt(gr.getSeqid(), gr.getStart(), gr.getEnd()).getBaseString();
                if (gr.getStrand().equals(GFF3Record.Strand.NEGATIVE)) {
                    seq = SequenceUtils.reverseComplement(seq);
                }

                out.println(">" + id);
                out.println(seq);
            }
        }
    }
}
