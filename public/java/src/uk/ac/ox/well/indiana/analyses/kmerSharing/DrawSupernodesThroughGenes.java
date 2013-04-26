package uk.ac.ox.well.indiana.analyses.kmerSharing;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import processing.pdf.PGraphicsPDF;
import uk.ac.ox.well.indiana.sketches.Sketch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.io.utils.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.awt.*;
import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class DrawSupernodesThroughGenes extends Sketch {
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
    private int geneHeight = 3;
    private int geneMargin = 5;
    private int geneLength = 0;
    private int legendMargin = 0;

    private Set<String> relatedSequences;
    private Map<String, String> genes;
    private Map<String, String> kmerToGene;
    private Map<String, Color> geneToColor;
    private Set<String> usedGenes;

    private Color[] generateColors(int n) {
        Color[] cols = new Color[n];

        for(int i = 0; i < n; i++) {
            cols[i] = Color.getHSBColor((float) i / (float) n, 0.85f, 1.0f);
        }

        return cols;
    }

    private void loadGenes() {
        genes = new TreeMap<String, String>();
        geneToColor = new HashMap<String, Color>();

        for (GFF3Record record : GFF) {
            if ("gene".equals(record.getType())) {
                String seq = new String(FASTA.getSubsequenceAt(record.getSeqid(), record.getStart(), record.getEnd()).getBases());

                String geneName = record.getAttribute("ID");

                genes.put(geneName, seq);
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
        usedGenes = new HashSet<String>();

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
                usedGenes.addAll(genesInSupernode);

                relatedSequences.add(supernode);
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

        for (String gene : usedGenes) {
            if (genes.get(gene).length() > geneLength) {
                geneLength = genes.get(gene).length();
            }
        }

        log.info("Used genes: {}", usedGenes);
    }

    public void setup() {
        int windowWidth = 2*horizontalMargin + geneLength;
        int windowHeight = 2*verticalMargin + geneHeight*(relatedSequences.size() + 3) + geneMargin*relatedSequences.size() + legendMargin;

        log.info("supernodes={}, width={} height={}", relatedSequences.size(), windowWidth, windowHeight);

        size(windowWidth, windowHeight, PGraphicsPDF.PDF, out.getAbsolutePath());

        for (String supernode : relatedSequences) {
            int sxpos0 = horizontalMargin;
            int sxpos1 = sxpos0 + supernode.length();
            int sypos = verticalMargin + geneHeight + geneMargin;

            stroke(Color.BLACK.getRGB());
            line(sxpos0, sypos, sxpos1, sypos);

            Map<String, Integer> genePositions = new HashMap<String, Integer>();
            Map<String, Integer> kmerPositions = new HashMap<String, Integer>();

            int i = 3;
            for (String gene : usedGenes) {
                String seq = genes.get(gene);

                int xpos0 = horizontalMargin;
                int xpos1 = xpos0 + seq.length();
                int ypos  = verticalMargin + (i*geneHeight) + (i*geneMargin);

                stroke(Color.BLACK.getRGB());
                line(xpos0, ypos, xpos1, ypos);

                genePositions.put(gene, ypos);
                for (int j = 0; j <= seq.length() - KMER_SIZE; j++) {
                    String fw = SequenceUtils.alphanumericallyLowestOrientation(seq.substring(j, j + KMER_SIZE));

                    kmerPositions.put(fw, xpos0 + j);
                }

                i++;
            }

            for (int j = 0; j <= supernode.length() - KMER_SIZE; j++) {
                String fw = SequenceUtils.alphanumericallyLowestOrientation(supernode.substring(j, j + KMER_SIZE));

                if (kmerToGene.containsKey(fw)) {
                    int kxpos = kmerPositions.get(fw);
                    int kypos = genePositions.get(kmerToGene.get(fw));

                    Color color = geneToColor.get(kmerToGene.get(fw));

                    stroke(color.getRGB());
                    line(kxpos, kypos - 1, kxpos, kypos + 1);
                }
            }

            PGraphicsPDF pdf = (PGraphicsPDF) g;
            pdf.nextPage();
        }

        exit();
    }
}
