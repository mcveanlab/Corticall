package uk.ac.ox.well.indiana.commands.prg;

import com.google.common.base.Joiner;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataFrame;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

@Description(text="Hierarchically cluster genes by sequence similarity")
public class clustergenes extends Module {
    @Argument(fullName="fasta", shortName="f", doc="ID:FASTA key-value pair")
    public HashMap<String, IndexedFastaSequenceFile> FASTAS;

    @Argument(fullName="gff", shortName="g", doc="ID:GFF key-value pair")
    public HashMap<String, GFF3> GFFS;

    @Argument(fullName="proteinKmerSize", shortName="pk", doc="Kmer size (in amino-acid space)")
    public Integer PROTEIN_KMER_SIZE = 7;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        // Make a map of all protein kmers in all genome
        log.info("Constructing map of all protein kmers in all genomes...");

        Map<String, Set<String>> kmerToSeqMap = new HashMap<String, Set<String>>();

        for (String id : FASTAS.keySet()) {
            IndexedFastaSequenceFile fasta = FASTAS.get(id);

            if (GFFS.containsKey(id)) {
                log.info("  Processing genome '{}'", id);

                GFF3 gff = GFFS.get(id);

                for (GFF3Record gr : gff) {
                    String geneName = gr.getAttribute("ID");

                    if ("gene".equals(gr.getType())) {
                        Collection<GFF3Record> exons = GFF3.getType("exon", gff.getChildren(gr));
                        String cds = SequenceUtils.extractCodingSequence(exons, fasta);

                        if (cds.contains("N")) {
                            log.warn("    Skipping gene {} (contains Ns in the sequence)", geneName);
                        } else {
                            String pseq = SequenceUtils.translateCodingSequence(cds);

                            for (int i = 0; i <= pseq.length() - PROTEIN_KMER_SIZE; i++) {
                                String kmer = pseq.substring(i, i + PROTEIN_KMER_SIZE);

                                if (!kmer.contains("N")) {
                                    if (!kmerToSeqMap.containsKey(kmer)) {
                                        kmerToSeqMap.put(kmer, new HashSet<String>());
                                    }

                                    kmerToSeqMap.get(kmer).add(geneName);
                                }
                            }
                        }
                    }
                }
            } else {
                log.warn("  Skipping genome '{}': did not find a GFF with the same ID, '{}'", fasta, id);
            }
        }

        // Compute sharing between all genes
        log.info("Computing sharing between all genes in all genomes...");

        DataFrame<String, String, Float> rmat = new DataFrame<String, String, Float>(0.0f);

        Map<String, Float> totalKmers = new HashMap<String, Float>();

        int kmersProcessed = 0;
        for (String kmer : kmerToSeqMap.keySet()) {
            if (kmersProcessed % (kmerToSeqMap.size() / 10) == 0) {
                log.info("  Processed {}/{} (~{}%) kmers", kmersProcessed, kmerToSeqMap.size(), String.format("%.1f", 100.0*kmersProcessed / kmerToSeqMap.size()));
            }
            kmersProcessed++;

            for (String name1 : kmerToSeqMap.get(kmer)) {
                for (String name2 : kmerToSeqMap.get(kmer)) {
                    rmat.set(name1, name2, rmat.get(name1, name2) + 1.0f);
                }

                if (!totalKmers.containsKey(name1)) {
                    totalKmers.put(name1, 0.0f);
                }
                totalKmers.put(name1, totalKmers.get(name1) + 1.0f);
            }
        }

        // Print table
        log.info("Writing table to disk...");

        Collection<String> rowNames = rmat.getRowNames();
        Collection<String> colNames = rmat.getColNames();

        out.println("\t" + Joiner.on("\t").join(colNames));

        for (String name1 : rowNames) {
            Collection<String> fields = new ArrayList<String>();
            fields.add(name1);

            for (String name2 : colNames) {
                float intersection = rmat.get(name1, name2);
                float total = intersection / (totalKmers.get(name1) < totalKmers.get(name2) ? totalKmers.get(name1) : totalKmers.get(name2));
                float distance = 1.0f - total;

                fields.add(String.valueOf(distance));
            }

            out.println(Joiner.on("\t").join(fields));
        }
    }
}
