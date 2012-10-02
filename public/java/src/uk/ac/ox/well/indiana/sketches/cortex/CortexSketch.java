package uk.ac.ox.well.indiana.sketches.cortex;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.sketches.Sketch;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexColor;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class CortexSketch extends Tool {
    @Argument(fullName="cortexGraph", shortName="cg", doc="A binary Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="genesFasta", shortName="gf", doc="A FASTA-format genes file")
    public FastaSequenceFile GENES_FASTA;

    private byte[] getReverseComplement(byte[] kmer) {
        byte[] rc = new byte[kmer.length];

        for (int i = 0; i < kmer.length; i++) {
            byte rcBase = 'N';
            switch (kmer[i]) {
                case 'A':
                    rcBase = 'T'; break;
                case 'C':
                    rcBase = 'G'; break;
                case 'G':
                    rcBase = 'C'; break;
                case 'T':
                    rcBase = 'A'; break;
            }

            rc[kmer.length - 1 - i] = rcBase;
        }

        return rc;
    }

    private byte[] getCortexCompatibleOrientation(byte[] kmer) {
        byte[] rc = getReverseComplement(kmer);

        String kmerStr = new String(kmer);
        String rcStr = new String(rc);

        return (kmerStr.compareTo(rcStr) < 0) ? kmer : rc;
    }

    private HashMap<String, String> loadGeneKmers(FastaSequenceFile genes, int kmerSize) {
        HashMap<String, String> kmerMap = new HashMap<String, String>();
        HashMap<String, Integer> kmerCoverageMap = new HashMap<String, Integer>();

        ReferenceSequence seq;
        while ((seq = genes.nextSequence()) != null) {
            String[] namePieces = seq.getName().split("\\s+");
            String geneName = namePieces[0];

            for (int i = 0; i < seq.length() - kmerSize; i++) {
                String kmer = new String(getCortexCompatibleOrientation(Arrays.copyOfRange(seq.getBases(), i, i + kmerSize)));

                kmerMap.put(kmer, geneName);

                if (!kmerCoverageMap.containsKey(kmer)) {
                    kmerCoverageMap.put(kmer, 1);
                } else {
                    kmerCoverageMap.put(kmer, kmerCoverageMap.get(kmer) + 1);
                }
            }
        }

        for (String kmer : kmerCoverageMap.keySet()) {
            if (kmerCoverageMap.get(kmer) > 1) {
                kmerMap.remove(kmer);
            }
        }

        return kmerMap;
    }

    public int execute() {
        HashMap<String, String> kmerMap = loadGeneKmers(GENES_FASTA, CORTEX_GRAPH.getKmerSize());

        System.out.println("Total kmers loaded: " + kmerMap.size());

        for (CortexRecord cr : CORTEX_GRAPH) {
            int color = CORTEX_GRAPH.getColorForSampleName("genes");

            if (color >= 0) {
                if (cr.getCoverages()[color] == 1) {
                    String kmer = cr.getKmerString();

                    if (kmerMap.containsKey(kmer)) {
                        String geneName = kmerMap.get(kmer);

                        System.out.println(geneName + " " + cr);
                    }
                }
            }
        }

        return 0;
    }

    /*
    public void setup() {
        HashMap<String, String> kmerMap = loadGeneKmers(GENES_FASTA, CORTEX_GRAPH.getKmerSize());

//        System.out.println(CORTEX_GRAPH);
//        System.out.println(GENES_FASTA);

//        for (CortexRecord cr : CORTEX_GRAPH) {
//            System.out.println(cr);
//        }

        size(400, 400);
    }

    public void draw() {
        if (mousePressed) {
            ellipse(mouseX, mouseY, 20, 20);
        }
    }
    */
}
