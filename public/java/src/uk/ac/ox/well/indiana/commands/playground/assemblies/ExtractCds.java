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
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
                //log.info("{}", gr);

                List<String> exons = new ArrayList<>();

                for (GFF3Record grcds : GFF3.getType("CDS", GFF.getContained(gr))) {
                    log.info("  {}", grcds);

                    ReferenceSequence rseq = REF.getSubsequenceAt(grcds.getSeqid(), grcds.getStart(), grcds.getEnd());
                    exons.add(rseq.getBaseString());
                }

                String cds = Joiner.on("").join(exons);
                if (gr.getStrand().equals(GFF3Record.Strand.NEGATIVE)) {
                    cds = SequenceUtils.reverseComplement(cds);
                }

                if (exons.size() > 0) {
                    out.println(">" + gr.getAttribute("ID") + ":" + (gr.getStrand().equals(GFF3Record.Strand.POSITIVE) ? "+" : "-"));
                    out.println(cds);
                    out.println(Joiner.on(" ").join(exons));
                }
            }
        }
    }
}
