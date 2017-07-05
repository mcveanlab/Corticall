package uk.ac.ox.well.indiana.commands.playground.assemblies;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by kiran on 05/07/2017.
 */
public class ExtractCds extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference")
    public IndexedFastaSequenceFile REF;

    @Argument(fullName="gff", shortName="g", doc="GFF")
    public GFF3 GFF;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        for (GFF3Record gr : GFF) {
            if (gr.getType().equals("gene")) {
                TreeSet<GFF3Record> grcds = new TreeSet<>(GFF3.getType("CDS", GFF.getContained(gr)));

                List<String> exons = new ArrayList<>();

                for (GFF3Record grcd : grcds) {
                    ReferenceSequence rseq = REF.getSubsequenceAt(grcd.getSeqid(), grcd.getStart(), grcd.getEnd());

                    if (gr.getStrand().equals(GFF3Record.Strand.NEGATIVE)) {
                        exons.add(0, SequenceUtils.reverseComplement(rseq.getBaseString()));
                    } else {
                        exons.add(rseq.getBaseString());
                    }
                }

                String cds = Joiner.on("").join(exons);

                if (exons.size() > 0) {
                    out.println(">" + gr.getAttribute("ID"));
                    out.println(cds);
                }
            }
        }
    }
}
