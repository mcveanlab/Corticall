package uk.ac.ox.well.indiana.tools.cortex;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.tools.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class CountGloballyUniqueGeneKmers extends Module {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="genesFasta", shortName="gf", doc="Genes fasta")
    public FastaSequenceFile GENES_FASTA;

    @Output
    public PrintStream out;

    private HashMap<String, String> loadGeneKmers(FastaSequenceFile genes, int kmerSize) {
        HashMap<String, String> kmerMap = new HashMap<String, String>();
        HashMap<String, Integer> kmerCoverageMap = new HashMap<String, Integer>();

        HashSet<String> geneNames = new HashSet<String>();

        ReferenceSequence seq;
        while ((seq = genes.nextSequence()) != null) {
            String[] namePieces = seq.getName().split("\\s+");
            String geneName = namePieces[0];

            geneNames.add(geneName);

            for (int i = 0; i < seq.length() - kmerSize; i++) {
                String kmer = new String(SequenceUtils.alphanumericallyLowestOrientation(Arrays.copyOfRange(seq.getBases(), i, i + kmerSize)));

                kmer = kmer.toUpperCase();
                if (!kmer.contains("N") && !kmer.contains(".")) {
                    kmerMap.put(kmer, geneName);

                    if (!kmerCoverageMap.containsKey(kmer)) {
                        kmerCoverageMap.put(kmer, 1);
                    } else {
                        kmerCoverageMap.put(kmer, kmerCoverageMap.get(kmer) + 1);
                    }
                }
            }
        }

        HashSet<String> genesInMap = new HashSet<String>();
        for (String kmer : kmerMap.keySet()) {
            genesInMap.add(kmerMap.get(kmer));
        }

        for (String kmer : kmerCoverageMap.keySet()) {
            if (kmerCoverageMap.get(kmer) > 1) {
                kmerMap.remove(kmer);
            }
        }

        genesInMap.clear();
        for (String kmer : kmerMap.keySet()) {
            genesInMap.add(kmerMap.get(kmer));
        }

        for (String geneName : geneNames) {
            if (!genesInMap.contains(geneName)) {
                System.out.println("Gene '" + geneName + "' was removed from the gene list.");
            }
        }

        return kmerMap;
    }

    @Override
    public void execute() {
        HashMap<String, String> kmerMap = loadGeneKmers(GENES_FASTA, CORTEX_GRAPH.getKmerSize());

        int colorGenes = CORTEX_GRAPH.getColorForSampleName("genes");
        int color3D7 = CORTEX_GRAPH.getColorForSampleName("3D7");

        HashMap<String, Integer> globallyUniqueKmers = new HashMap<String, Integer>();
        HashMap<String, Integer> totalKmers = new HashMap<String, Integer>();

        for (CortexRecord cr : CORTEX_GRAPH) {
            if (cr.getCoverages()[colorGenes] == 1 && cr.getCoverages()[color3D7] == 1) {
                if (kmerMap.containsKey(cr.getKmerAsString())) {
                    String geneName = kmerMap.get(cr.getKmerAsString());

                    boolean globallyUnique = true;

                    for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
                        if (cr.getCoverages()[color] != 1) {
                            globallyUnique = false;
                            break;
                        }
                    }

                    if (!totalKmers.containsKey(geneName)) {
                        totalKmers.put(geneName, 1);
                        globallyUniqueKmers.put(geneName, 0);
                    } else {
                        totalKmers.put(geneName, totalKmers.get(geneName) + 1);
                    }

                    if (globallyUnique) {
                        globallyUniqueKmers.put(geneName, globallyUniqueKmers.get(geneName) + 1);
                    }
                }
            }
        }

        out.printf("gene totalKmers.%d globallyUniqueKmers.%d\n", CORTEX_GRAPH.getKmerSize(), CORTEX_GRAPH.getKmerSize());
        for (String geneName : totalKmers.keySet()) {
            out.printf("%s %d %d\n", geneName, totalKmers.get(geneName), globallyUniqueKmers.get(geneName));
        }
    }
}
