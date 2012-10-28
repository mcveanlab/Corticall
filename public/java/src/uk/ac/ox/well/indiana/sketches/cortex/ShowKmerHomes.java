package uk.ac.ox.well.indiana.sketches.cortex;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import processing.core.PFont;
import uk.ac.ox.well.indiana.sketches.Sketch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.util.HashMap;
import java.util.TreeMap;

public class ShowKmerHomes extends Sketch {
    @Argument(fullName="fasta", shortName="f", doc="Reference fasta file")
    public FastaSequenceFile REFERENCE;

    @Argument(fullName="gene", shortName="g", doc="The gene fasta to examine")
    public FastaSequenceFile GENES;

    public HashMap<String, String> loadGenes(FastaSequenceFile genesFasta) {
        HashMap<String, String> genes = new HashMap<String, String>();

        ReferenceSequence gene;
        while ((gene = genesFasta.nextSequence()) != null) {
            String[] names = gene.getName().split("\\s+");
            genes.put(new String(gene.getBases()), names[0]);
        }

        return genes;
    }

    public void setup() {
        size(600, 600);
        //noLoop();

        background(255);

        PFont f = createFont("Helvetica", 16, true);
        textFont(f, 16);
        fill(0);

        int maxSeqLength = 0;

        HashMap<String, String> geneSeqs = loadGenes(GENES);
        TreeMap<String, ReferenceSequence> contigs = new TreeMap<String, ReferenceSequence>();

        ReferenceSequence seq;
        while ((seq = REFERENCE.nextSequence()) != null) {
            if (seq.length() > maxSeqLength) {
                maxSeqLength = seq.length();
            }

            String[] names = seq.getName().split("\\s+");

            contigs.put(names[0], seq);
        }

        int margin = 10;
        float verticalStep = ((float) (height - 2*margin))/contigs.size();

        int i = 0;
        for (String contigName : contigs.keySet()) {
            ReferenceSequence contig = contigs.get(contigName);
            String contigSeq = new String(contig.getBases());

            float y = 2*margin + (i*verticalStep);
            float x2 = margin + contig.length() * (width - 2*margin) / maxSeqLength;

            line(margin, y, x2, y);

            text(contigName, margin, y - 2);

            for (String geneSeq : geneSeqs.keySet()) {
                System.out.println(geneSeq);

                String rc = new String(SequenceUtils.getReverseComplement(geneSeq.getBytes()));
                if (contigSeq.contains(rc)) {
                    int x3 = margin + contigSeq.indexOf(rc) * (width - 2*margin) / maxSeqLength;
                    int x4 = margin + (contigSeq.indexOf(rc) + geneSeq.length()) * (width - 2*margin) / maxSeqLength;

                    line(x3, y + 2, x4, y);

                    System.out.println(x3 + " " + x4);
                }
            }

            i++;
        }
    }
}
