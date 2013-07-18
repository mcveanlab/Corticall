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
import java.util.TreeMap;

public class CountGeneKmerCoverageByPosition extends Module {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="genesFasta", shortName="gf", doc="Genes fasta file")
    public FastaSequenceFile GENES_FASTA;

    @Argument(fullName="geneName", shortName="gn", doc="Gene name")
    public String GENE_NAME;

    @Output
    public PrintStream out;

    private class GeneNameAndPosition {
        public String geneName;
        public int pos;
    }

    private HashMap<Integer, GeneNameAndPosition> loadGeneKmers(FastaSequenceFile genesFasta, int kmerSize) {
        HashMap<Integer, GeneNameAndPosition> geneKmers = new HashMap<Integer, GeneNameAndPosition>();

        ReferenceSequence seq;
        while ((seq = genesFasta.nextSequence()) != null) {
            String[] names = seq.getName().split("\\s+");
            String name = names[0];

            if (name.contains(GENE_NAME)) {
                byte[] b = seq.getBases();

                for (int i = 0; i < b.length - kmerSize; i++) {
                    String kmer = new String(SequenceUtils.alphanumericallyLowestOrientation(Arrays.copyOfRange(b, i, i + kmerSize)));

                    GeneNameAndPosition gnp = new GeneNameAndPosition();
                    gnp.geneName = name;
                    gnp.pos = i;

                    geneKmers.put(kmer.hashCode(), gnp);
                }
            }
        }

        return geneKmers;
    }

    @Override
    public void execute() {
        HashMap<Integer, GeneNameAndPosition> geneKmers = loadGeneKmers(GENES_FASTA, CORTEX_GRAPH.getKmerSize());
        HashMap<String, TreeMap<Integer, HashMap<String, Integer>>> kmerCoverage = new HashMap<String, TreeMap<Integer, HashMap<String, Integer>>>();

        int colorGenes = CORTEX_GRAPH.getColorForSampleName("genes");
        int color3D7 = CORTEX_GRAPH.getColorForSampleName("3D7");

        if (colorGenes >= 0 && color3D7 >= 0) {
            int numRecords = 0;
            for (CortexRecord cr : CORTEX_GRAPH) {
                if (numRecords % (int) (CORTEX_GRAPH.getNumRecords()/10) == 0) {
                    System.out.println("Saw " + numRecords + "/" + CORTEX_GRAPH.getNumRecords() + " records");
                }
                numRecords++;

                if (cr.getCoverages()[colorGenes] == 1 && cr.getCoverages()[color3D7] == 1) {
                    String kmer = cr.getKmerAsString();

                    if (geneKmers.containsKey(kmer.hashCode())) {
                        GeneNameAndPosition gnp = geneKmers.get(kmer.hashCode());

                        if (!kmerCoverage.containsKey(gnp.geneName)) {
                            kmerCoverage.put(gnp.geneName, new TreeMap<Integer, HashMap<String, Integer>>());
                        }

                        if (!kmerCoverage.get(gnp.geneName).containsKey(gnp.pos)) {
                            kmerCoverage.get(gnp.geneName).put(gnp.pos, new HashMap<String, Integer>());
                        }

                        for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
                            String sampleName = CORTEX_GRAPH.getColor(color).getSampleName();

                            kmerCoverage.get(gnp.geneName).get(gnp.pos).put(sampleName, cr.getCoverages()[color]);
                        }
                    }
                }
            }
        }

        for (String geneName : kmerCoverage.keySet()) {
            TreeMap<Integer, HashMap<String, Integer>> positionCoverage = kmerCoverage.get(geneName);

            String header = "gene sample";
            for (int pos : positionCoverage.keySet()) {
                header += " " + pos;
            }

            out.println(header);

            for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
                String sampleName = CORTEX_GRAPH.getColor(color).getSampleName();

                String line = geneName + " " + sampleName;

                for (int pos : positionCoverage.keySet()) {
                    line += " " + positionCoverage.get(pos).get(sampleName);
                }

                out.println(line);
            }
        }
    }
}
