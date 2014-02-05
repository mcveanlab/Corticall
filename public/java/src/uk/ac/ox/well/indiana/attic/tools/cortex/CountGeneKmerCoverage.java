package uk.ac.ox.well.indiana.attic.tools.cortex;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexColor;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;

public class CountGeneKmerCoverage extends Module {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="genesFasta", shortName="gf", doc="Genes fasta file")
    public FastaSequenceFile GENES_FASTA;

    @Output
    public PrintStream out;

    private HashMap<Integer, String> loadGeneKmers(FastaSequenceFile genesFasta, int kmerSize) {
        HashMap<Integer, String> geneKmers = new HashMap<Integer, String>();

        ReferenceSequence seq;
        while ((seq = genesFasta.nextSequence()) != null) {
            String[] names = seq.getName().split("\\s+");
            String name = names[0];

            byte[] b = seq.getBases();

            for (int i = 0; i < b.length - kmerSize; i++) {
                String kmer = new String(SequenceUtils.alphanumericallyLowestOrientation(Arrays.copyOfRange(b, i, i + kmerSize)));

                geneKmers.put(kmer.hashCode(), name);
            }
        }

        return geneKmers;
    }

    private class GeneCoverage {
        public int numKmers = 0;

        public HashMap<String, Integer> cov0;
        public HashMap<String, Integer> cov1;
        public HashMap<String, Integer> cov2;
        public HashMap<String, Integer> covMore;

        public GeneCoverage() {
            cov0 = new HashMap<String, Integer>();
            cov1 = new HashMap<String, Integer>();
            cov2 = new HashMap<String, Integer>();
            covMore = new HashMap<String, Integer>();

            for (CortexColor cc : CORTEX_GRAPH.getColors()) {
                cov0.put(cc.getSampleName(), 0);
                cov1.put(cc.getSampleName(), 0);
                cov2.put(cc.getSampleName(), 0);
                covMore.put(cc.getSampleName(), 0);
            }
        }
    }

    @Override
    public void execute() {
        HashMap<Integer, String> geneKmers = loadGeneKmers(GENES_FASTA, CORTEX_GRAPH.getKmerSize());

        HashMap<String, GeneCoverage> geneCoverage = new HashMap<String, GeneCoverage>();

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
                        String geneName = geneKmers.get(kmer.hashCode());

                        if (!geneCoverage.containsKey(geneName)) {
                            geneCoverage.put(geneName, new GeneCoverage());
                        }

                        geneCoverage.get(geneName).numKmers++;

                        for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
                            String sampleName = CORTEX_GRAPH.getColor(color).getSampleName();
                            int coverage = cr.getCoverages()[color];

                            if (coverage == 0) {
                                int oldCount = geneCoverage.get(geneName).cov0.get(sampleName);
                                geneCoverage.get(geneName).cov0.put(sampleName, oldCount + 1);
                            } else if (coverage == 1) {
                                int oldCount = geneCoverage.get(geneName).cov1.get(sampleName);
                                geneCoverage.get(geneName).cov1.put(sampleName, oldCount + 1);
                            } else if (coverage == 2) {
                                int oldCount = geneCoverage.get(geneName).cov2.get(sampleName);
                                geneCoverage.get(geneName).cov2.put(sampleName, oldCount + 1);
                            } else {
                                int oldCount = geneCoverage.get(geneName).covMore.get(sampleName);
                                geneCoverage.get(geneName).covMore.put(sampleName, oldCount + 1);
                            }
                        }
                    }
                }
            }
        } else {
            // yell at user
        }

        String header = "gene kmerSize numKmers coverage";
        for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
            String sampleName = CORTEX_GRAPH.getColor(color).getSampleName();

            header += " " + sampleName;
        }

        out.println(header);

        for (String geneName : geneCoverage.keySet()) {
            GeneCoverage gc = geneCoverage.get(geneName);

            for (int coverage = 0; coverage <= 3; coverage++) {
                String value = geneName + " " + CORTEX_GRAPH.getKmerSize() + " " + gc.numKmers + " " + coverage;

                HashMap<String, Integer> coverageCounts;
                switch (coverage) {
                    case 0:
                        coverageCounts = gc.cov0; break;
                    case 1:
                        coverageCounts = gc.cov1; break;
                    case 2:
                        coverageCounts = gc.cov2; break;
                    default:
                        coverageCounts = gc.covMore; break;
                }

                for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
                    String sampleName = CORTEX_GRAPH.getColor(color).getSampleName();

                    value += " " + coverageCounts.get(sampleName);
                }

                out.println(value);
            }
        }
    }
}
