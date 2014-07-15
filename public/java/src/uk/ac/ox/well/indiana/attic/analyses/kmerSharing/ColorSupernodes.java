package uk.ac.ox.well.indiana.attic.analyses.kmerSharing;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import processing.pdf.PGraphicsPDF;
import uk.ac.ox.well.indiana.commands.Sketch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.awt.*;
import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class ColorSupernodes extends Sketch {
    @Argument(fullName="reference", shortName="R", doc="Reference genome")
    public IndexedFastaSequenceFile FASTA;

    @Argument(fullName="gff", shortName="gff", doc="GFF of regions of interest")
    public GFF3 GFF;

    @Argument(fullName="relatedSequences", shortName="rel", doc="Related sequences file")
    public File RELATED_SEQUENCES;

    @Argument(fullName="kmerReferencePanel", shortName="krp", doc="Kmer reference panel")
    public File KMER_REFERENCE_PANEL;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Argument(fullName="geneOfInterest", shortName="goi", doc="Genes of interest")
    public String GENE_OF_INTEREST = "PF3D7_0100100";

    @Output(fullName="fastaOut", shortName="fo", doc="Fasta out")
    public PrintStream fout;

    @Output
    public File out;

    private int verticalMargin = 10;
    private int horizontalMargin = 10;
    private int supernodeHeight = 3;
    private int supernodeLength = 0;
    private int legendMargin = 60*supernodeHeight;

    private Set<String> relatedSequences;
    private Map<String, String> genes;
    private Map<String, String> kmerToGene;
    private Map<String, Color> geneToColor;

    private Color[] generateColors(int n) {
        Color[] cols = new Color[n];

        for(int i = 0; i < n; i++) {
            cols[i] = Color.getHSBColor((float) i / (float) n, 0.85f, 1.0f);
        }

        return cols;
    }

    private void loadGenes() {
        genes = new HashMap<String, String>();
        geneToColor = new HashMap<String, Color>();

        for (GFF3Record record : GFF) {
            if ("gene".equals(record.getType())) {
                String seq = new String(FASTA.getSubsequenceAt(record.getSeqid(), record.getStart(), record.getEnd()).getBases());

                String geneName = record.getAttribute("ID");

                genes.put(geneName, seq);

                /*
                for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                    String fw = SequenceUtils.alphanumericallyLowestOrientation(seq.substring(i, i + KMER_SIZE));

                    if (!kmerToGene.containsKey(fw)) {
                        kmerToGene.put(fw, geneName);
                    }
                }
                */
            }
        }

        Color[] colors = generateColors(genes.size());
        int i = 0;
        for (String gene : genes.keySet()) {
            geneToColor.put(gene, colors[i]);

            i++;
        }
    }

    private void loadRelatedSequences() {
        relatedSequences = new HashSet<String>();

        TableReader tr = new TableReader(RELATED_SEQUENCES);

        for (Map<String, String> te : tr) {
            String supernode = SequenceUtils.alphanumericallyLowestOrientation(te.get("superNode"));

            Set<String> genesInSupernode = new HashSet<String>();

            for (int i = 0; i <= supernode.length() - KMER_SIZE; i++) {
                String fw = SequenceUtils.alphanumericallyLowestOrientation(supernode.substring(i, i + KMER_SIZE));

                if (kmerToGene.containsKey(fw)) {
                    genesInSupernode.add(kmerToGene.get(fw));
                }
            }

            if (genesInSupernode.size() > 1 && genesInSupernode.contains(GENE_OF_INTEREST)) {
                log.debug("genesInSupernode: {}", genesInSupernode);

                relatedSequences.add(supernode);

                int length = supernode.length();

                if (length > supernodeLength) {
                    supernodeLength = length;
                }
            }
        }
    }

    private void loadKmerReferencePanel() {
        kmerToGene = new HashMap<String, String>();

        TableReader tr = new TableReader(KMER_REFERENCE_PANEL);

        for (Map<String, String> te : tr) {
            String kmer = te.get("kmer");
            String gene = te.get("genes");

            kmerToGene.put(kmer, gene);
        }
    }

    public void initialize() {
        loadGenes();
        loadKmerReferencePanel();
        loadRelatedSequences();
    }

    public void setup() {
        int windowWidth = 2*horizontalMargin + supernodeLength;
        int windowHeight = 2*verticalMargin + supernodeHeight*relatedSequences.size() + legendMargin;

        log.info("supernodes={}, width={} height={}", relatedSequences.size(), windowWidth, windowHeight);

        size(windowWidth, windowHeight, PGraphicsPDF.PDF, out.getAbsolutePath());

        Set<String> usedGenes = new HashSet<String>();

        int i = 0;
        for (String supernode : relatedSequences) {
            String supernodeRc = SequenceUtils.reverseComplement(supernode);

            fout.printf(">sn.%d_sn.%d\n", supernode.hashCode(), supernodeRc.hashCode());
            fout.printf("%s\n", supernode);

            int xpos0 = horizontalMargin;
            int xpos1 = xpos0 + supernode.length();
            int ypos  = verticalMargin + (i*supernodeHeight);

            stroke(Color.BLACK.getRGB());
            line(xpos0, ypos, xpos1, ypos);

            i++;

            for (int j = 0; j <= supernode.length() - KMER_SIZE; j++) {
                String fw = SequenceUtils.alphanumericallyLowestOrientation(supernode.substring(j, j + KMER_SIZE));

                if (kmerToGene.containsKey(fw)) {
                    Color color = geneToColor.get(kmerToGene.get(fw));

                    int kxpos = xpos0 + j;
                    int kypos0 = ypos - 1;
                    int kypos1 = ypos + 1;

                    stroke(color.getRGB());
                    line(kxpos, kypos0, kxpos, kypos1);

                    usedGenes.add(kmerToGene.get(fw));
                }
            }
        }

        i = 0;
        for (String gene : usedGenes) {
            int xpos = horizontalMargin;
            int ypos = verticalMargin + supernodeHeight*relatedSequences.size() + (i*supernodeHeight);
            int ypos0 = ypos - 1;
            int ypos1 = ypos + 1;

            Color color = geneToColor.get(gene);

            stroke(color.getRGB());
            line(xpos, ypos0, xpos, ypos1);

            textSize(2);
            fill(Color.BLACK.getRGB());
            text(gene, xpos + 2, ypos);

            i++;
        }

        exit();
    }
}
