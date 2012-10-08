package uk.ac.ox.well.indiana.tools.cortex;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashSet;

public class PrintGeneKmersInCortexFile extends Tool {
    @Argument(fullName="cortexGraph", shortName="cg", doc="The Cortex binary to read")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="genesFasta", shortName="gf", doc="The genes fasta file")
    public FastaSequenceFile GENES_FASTA;

    @Argument(fullName="gene", shortName="g", doc="The gene to examine")
    public String GENE_NAME;

    @Output
    public PrintStream out;

    private HashSet<String> getGeneKmers() {
        HashSet<String> geneKmers = new HashSet<String>();

        ReferenceSequence seq;
        while ((seq = GENES_FASTA.nextSequence()) != null) {
            String[] pieces = seq.getName().split("\\s+");
            String geneName = pieces[0];

            if (geneName.equalsIgnoreCase(GENE_NAME)) {
                for (int i = 0; i < seq.length() - CORTEX_GRAPH.getKmerSize(); i++) {
                    String kmer = new String(SequenceUtils.getCortexCompatibleOrientation(Arrays.copyOfRange(seq.getBases(), i, i + CORTEX_GRAPH.getKmerSize())));

                    geneKmers.add(kmer);
                }
            }
        }

        return geneKmers;
    }

    @Override
    public int execute() {
        HashSet<String> geneKmers = getGeneKmers();
        HashSet<String> seenKmers = new HashSet<String>();

        int i = 0;
        for (CortexRecord cr : CORTEX_GRAPH) {
            String kmer = cr.getKmerString();

            if (geneKmers.contains(kmer)) {
                seenKmers.add(kmer);

                //out.println(GENE_NAME + " " + cr);
            }

            i++;
            if (i % 100000 == 0) {
                System.out.println("Read " + i + " Cortex records");
            }
        }

        out.println("     gene kmers: " + geneKmers.size());
        out.println("    graph kmers: " + seenKmers.size());

        geneKmers.removeAll(seenKmers);

        out.println("remaining kmers: " + geneKmers.size());

        for (String kmer : geneKmers) {
            out.println("\t" + kmer);
        }

        return 0;
    }
}
