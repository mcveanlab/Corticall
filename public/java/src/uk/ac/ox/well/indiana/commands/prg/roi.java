package uk.ac.ox.well.indiana.commands.prg;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;

@Description(text="Extracts regions of interest from a FASTA file")
public class roi extends Module {
    @Argument(fullName="fasta", shortName="f", doc="ID:FASTA key-value pair")
    public HashMap<String, IndexedFastaSequenceFile> FASTAS;

    @Argument(fullName="gff", shortName="g", doc="ID:GFF key-value pair")
    public HashMap<String, GFF3> GFFS;

    @Argument(fullName="ids", shortName="i", doc="Gene IDs to extract from GFF file", required=false)
    public HashSet<String> IDS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        log.info("Extracting {} genes of interest from FASTA files...", IDS.size());

        int seqsTotal = 0;
        for (String id : FASTAS.keySet()) {
            IndexedFastaSequenceFile fasta = FASTAS.get(id);

            if (GFFS.containsKey(id)) {
                log.info("  Processing genome '{}'", id);

                GFF3 gff = GFFS.get(id);

                int seqsFound = 0;
                for (GFF3Record gr : gff) {
                    if ("gene".equalsIgnoreCase(gr.getType())) {
                        String gid = gr.getAttribute("ID");

                        if (IDS == null || IDS.contains(gid)) {
                            String gene = SequenceUtils.extractGeneSequence(gr, fasta);

                            if (gene.contains("N")) {
                                log.warn("    Skipping gene {} (contains Ns in the sequence)", gid);
                            } else {
                                seqsFound++;

                                out.println(">" + gid + ".gene\n" + gene);

                                Collection<GFF3Record> exons = GFF3.getType("exon", gff.getChildren(gr));
                                if (!exons.isEmpty()) {
                                    String cds = SequenceUtils.extractCodingSequence(exons, fasta);
                                    String tr  = SequenceUtils.translateCodingSequence(cds);

                                    out.println(">" + gid + ".cds\n" + cds);
                                    out.println(">" + gid + ".tr\n" + tr);
                                }
                            }
                        }
                    }
                }

                log.info("    Extracted {} genes", seqsFound);
                seqsTotal += seqsFound;
            }
        }

        log.info("Extracted {} genes total", seqsTotal);
    }
}
