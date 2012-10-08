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
import java.util.*;

public class SortCortexRecords extends Tool {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="genesSampleName", shortName="gsn", doc="Sample name for genes color in Cortex file")
    public String GENES_SAMPLE_NAME = "genes";

    @Argument(fullName="genesFasta", shortName="gf", doc="Genes fasta file")
    public FastaSequenceFile GENES_FASTA;

    @Argument(fullName="gene", shortName="g", doc="Gene to process")
    public String GENE_NAME = null;

    @Output
    public PrintStream out;

    private HashMap<String, String> loadGenesMap(FastaSequenceFile fasta, String geneName, int kmerSize) {
        HashMap<String, String> geneMap = new HashMap<String, String>();

        ReferenceSequence seq;
        while ((seq = fasta.nextSequence()) != null) {
            String[] names = seq.getName().split("\\s+");

            if (GENE_NAME == null || names[0].equalsIgnoreCase(GENE_NAME)) {
                for (int i = 0; i < seq.length() - kmerSize; i++) {
                    String kmer = new String(SequenceUtils.getCortexCompatibleOrientation(Arrays.copyOfRange(seq.getBases(), i, i+kmerSize)));

                    geneMap.put(kmer, names[0]);
                }
            }
        }

        return geneMap;
    }

    @Override
    public int execute() {
        HashMap<String, String> geneMap = loadGenesMap(GENES_FASTA, GENE_NAME, CORTEX_GRAPH.getKmerSize());

        TreeSet<CortexRecord> sortedRecords = new TreeSet<CortexRecord>();

        long recordsSeen = 0;
        long recordsWritten = 0;
        long recordsTotal = CORTEX_GRAPH.getNumRecords();

        int genesColor = CORTEX_GRAPH.getColorForSampleName(GENES_SAMPLE_NAME);

        for (CortexRecord cr : CORTEX_GRAPH) {
            if (cr.getCoverages()[genesColor] > 0) {
                if (geneMap.containsKey(cr.getKmerString())) {
                    sortedRecords.add(cr);
                }

                recordsSeen++;
                if (recordsSeen % 100000 == 0) {
                    System.out.println("records seen: " + recordsSeen + "/" + recordsTotal);
                }
            }
        }

        System.out.println("loaded " + recordsSeen + "/" + recordsTotal);

        for (CortexRecord cr : sortedRecords) {
            out.println(cr);

            recordsWritten++;
            if (recordsWritten % 100000 == 0) {
                System.out.println("records written: " + recordsWritten + "/" + recordsTotal);
            }
        }

        System.out.println("wrote " + recordsWritten + "/" + recordsTotal);

        return 0;
    }
}
