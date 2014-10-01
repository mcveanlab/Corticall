package uk.ac.ox.well.indiana.attic.analyses.PacBio;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;

import java.io.PrintStream;

public class ExtractSequenceFeature extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference file (FASTA)")
    public IndexedFastaSequenceFile FASTA;

    @Argument(fullName="gff", shortName="g", doc="GFF file")
    public GFF3 GFF;

    @Argument(fullName="type", shortName="t", doc="Type to select(gene or exon)")
    public String TYPE;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        for (GFF3Record gr : GFF) {
            if (gr.getType().equals(TYPE)) {
                ReferenceSequence rseq = FASTA.getSubsequenceAt(gr.getSeqid(), gr.getStart(), gr.getEnd());
                String seq = new String(rseq.getBases());
                String id = gr.getAttribute("ID");

                out.println(">" + id);
                out.println(seq);
            }
        }
    }
}
