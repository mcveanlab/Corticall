package uk.ac.ox.well.indiana.commands.playground.geneannotation;

import com.google.api.services.genomics.model.Reference;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.BwaAligner;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.LastzAligner;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class ExtractExonSequences extends Module {
    @Argument(fullName="fasta", shortName="f", doc="Source FASTA file")
    public IndexedFastaSequenceFile SOURCE_FASTA;

    @Argument(fullName="gff", shortName="g", doc="Source GFF file")
    public GFF3 SOURCE_GFF;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        log.info("Loading gene and exon sequences...");

        int numGenes = 0;
        int numExons = 0;

        for (GFF3Record gr : SOURCE_GFF) {
            if (gr.getType().equals("gene")) {
                Collection<GFF3Record> exons = SOURCE_GFF.getChildren(gr);

                for (GFF3Record exon : exons) {
                    ReferenceSequence seq = new ReferenceSequence(
                            "genename_" + gr.getAttribute("ID") + "_exonname_" + exon.getAttribute("ID"),
                            numExons++,
                            SOURCE_FASTA.getSubsequenceAt(exon.getSeqid(), exon.getStart(), exon.getEnd()).getBases()
                    );

                    out.println(">" + seq.getName());
                    out.println(seq.getBaseString());
                }

                numGenes++;
            }
        }

        log.info("  {} genes, {} exons", numGenes, numExons);
    }
}
