package uk.ac.ox.well.indiana.attic.analyses.reconstruction;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import processing.pdf.PGraphicsPDF;
import uk.ac.ox.well.indiana.commands.Sketch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.awt.*;
import java.io.File;
import java.util.*;
import java.util.List;

public class VisualizeMultiPanelContigs extends Sketch {
    @Argument(fullName="contigTable", shortName="ct", doc="Contig table")
    public File CONTIG_TABLE;

    @Argument(fullName="geneOrder", shortName="go", doc="Gene order")
    public ArrayList<String> GENE_ORDER;

    @Argument(fullName="onlyTheseGenes", shortName="g", doc="Process only these genes", required=false)
    public HashSet<String> ONLY_THESE_GENES;

    @Argument(fullName="onlyTheseSamples", shortName="s", doc="Process only these samples", required=false)
    public HashSet<String> ONLY_THESE_SAMPLES;

    @Argument(fullName="gff", shortName="gff", doc="GFF file")
    public GFF3 GFF;

    @Argument(fullName="reference", shortName="R", doc="Reference FASTA file")
    public IndexedFastaSequenceFile FASTA;

    @Output
    public File out;

    //          sample      gene    entry
    private Map<String, Map<String, List<Map<String, String>>>> data;

    //          gene    color
    private Map<String, Color> geneColors;

    //          gene        kmers
    private Map<String, Set<CortexKmer>> geneKmers = new HashMap<String, Set<CortexKmer>>();

    private final int marginXRight = 50;
    private final int marginLabel = 100;
    private final int marginX = 10;
    private final int marginY = 10;
    private final int marginTitle = 30;
    private final int marginContig = 5;
    private final int contigHeight = 10;
    private int longestContig = 0;
    private int maxContigs = 0;

    public void initialize() {
        log.info("Process only these genes: {}", ONLY_THESE_GENES);
        log.info("Process only these samples: {}", ONLY_THESE_SAMPLES);

        data = new HashMap<String, Map<String, List<Map<String, String>>>>();
        geneColors = new TreeMap<String, Color>();

        TableReader tr = new TableReader(CONTIG_TABLE);

        for (Map<String, String> te : tr) {
            String sample = te.get("sample");
            Set<String> genes = new TreeSet<String>(Arrays.asList(te.get("genes").split(",")));

            if (genes.size() > 1) {
                if (!data.containsKey(sample)) {
                    data.put(sample, new TreeMap<String, List<Map<String, String>>>());
                }

                for (String gene : genes) {
                    if (!data.get(sample).containsKey(gene)) {
                        data.get(sample).put(gene, new ArrayList<Map<String, String>>());
                    }

                    data.get(sample).get(gene).add(te);
                    geneColors.put(gene, Color.BLACK);
                }
            }
        }

        int kmerSize = 0;
        for (String sample : data.keySet()) {
            for (String gene : data.get(sample).keySet()) {
                if ( (ONLY_THESE_SAMPLES == null || ONLY_THESE_SAMPLES.contains(sample)) && (ONLY_THESE_GENES == null || ONLY_THESE_GENES.contains(gene)) ) {
                    if (data.get(sample).get(gene).size() > maxContigs) {
                        maxContigs = data.get(sample).get(gene).size();
                    }

                    for (Map<String, String> te : data.get(sample).get(gene)) {
                        String contig = te.get("contig");
                        if (kmerSize == 0) {
                            String[] kmers = te.get("kmers").split(",");

                            kmerSize = kmers[0].length();
                        }

                        if (contig.length() > longestContig) {
                            longestContig = contig.length();
                        }
                    }
                }
            }
        }

        Color[] colors = generateColors(GENE_ORDER.size());

        int colorIndex = 0;
        for (String gene : GENE_ORDER) {
            geneColors.put(gene, colors[colorIndex]);

            colorIndex++;
        }

        for (String gene : ONLY_THESE_GENES) {
            GFF3Record record = GFF.getRecord(gene);

            String seq = new String(FASTA.getSubsequenceAt(record.getSeqid(), record.getStart(), record.getEnd()).getBases());
            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                CortexKmer kmer = new CortexKmer(seq.substring(i, i + kmerSize));

                if (!geneKmers.containsKey(gene)) {
                    geneKmers.put(gene, new HashSet<CortexKmer>());
                }

                geneKmers.get(gene).add(kmer);
            }
        }
    }

    // Taken from http://stackoverflow.com/questions/223971/how-to-generate-spectrum-color-palettes/223981#223981
    private Color[] generateColors(int n) {
        Color[] cols = new Color[n];

        for(int i = 0; i < n; i++) {
            cols[i] = Color.getHSBColor((float) i / (float) n, 0.85f, 1.0f);
        }

        return cols;
    }

    public void setup() {
        size(marginLabel + 2*marginX + longestContig + marginXRight, marginTitle + 2*marginY + (marginContig + contigHeight)*maxContigs, PGraphicsPDF.PDF, out.getAbsolutePath());

        background(Color.WHITE.getRGB());

        PGraphicsPDF pdf = (PGraphicsPDF) g;

        log.info("{}x{}", getWidth(), getHeight());

        int index = 0;
        for (String sample : data.keySet()) {
            for (String gene : data.get(sample).keySet()) {
                if ( (ONLY_THESE_SAMPLES == null || ONLY_THESE_SAMPLES.contains(sample)) && (ONLY_THESE_GENES == null || ONLY_THESE_GENES.contains(gene)) ) {
                    if (index > 0) { pdf.nextPage(); }
                    index++;

                    log.info("{} {} {}", sample, gene, data.get(sample).get(gene).size());

                    //background(204, 204, 204);
                    background(Color.WHITE.getRGB());

                    stroke(geneColors.get(gene).getRGB());
                    fill(geneColors.get(gene).getRGB(), 150.0f);
                    rect(marginX, marginTitle/2, 20, 7);

                    fill(Color.BLACK.getRGB());
                    textAlign(LEFT, CENTER);
                    textSize(15);
                    text("gene=" + gene + " sample=" + sample, marginX + 30, marginTitle/2);

                    for (int i = 0; i < data.get(sample).get(gene).size(); i++) {
                        Map<String, String> te = data.get(sample).get(gene).get(i);
                        String contig = te.get("contig");

                        Map<CortexKmer, String> kmersAndGenes = new HashMap<CortexKmer, String>();
                        String[] kmers = te.get("kmers").split(",");
                        String[] genes = te.get("genes").split(",");

                        int kmerSize = 0;
                        for (int j = 0; j < kmers.length; j++) {
                            kmersAndGenes.put(new CortexKmer(kmers[j]), genes[j]);

                            if (kmerSize == 0) {
                                kmerSize = kmers[j].length();
                            }
                        }

                        int xpos = marginLabel + marginX;
                        int ypos = marginTitle + marginY + i*(marginContig + contigHeight);

                        stroke(Color.BLACK.getRGB());
                        fill(Color.WHITE.getRGB());
                        rect(xpos - 1, ypos, contig.length() + 1, contigHeight);

                        fill(Color.BLACK.getRGB());
                        textAlign(LEFT, TOP);
                        textSize(10);
                        text(contig.length() + " bp", xpos + contig.length() + 3, ypos - 2);

                        boolean contigContainsKmerFromGene = false;

                        for (int k = 0; k <= contig.length() - kmerSize; k++) {
                            CortexKmer kmer = new CortexKmer(contig.substring(k, k+kmerSize));

                            if (kmersAndGenes.containsKey(kmer)) {
                                String geneName = kmersAndGenes.get(kmer);

                                if (!gene.equalsIgnoreCase(geneName) && ONLY_THESE_GENES.contains(geneName)) {
                                    contigContainsKmerFromGene = true;
                                }

                                byte[] kmerBytes = kmer.getKmerAsBytes();

                                int kmerHeight = contigHeight;

                                for (byte[] e1 : SequenceUtils.generateSequencesWithEditDistance1(kmerBytes)) {
                                    CortexKmer newCortexKmer = new CortexKmer(e1);

                                    for (String altGene : geneKmers.keySet()) {
                                        if (geneKmers.get(altGene).contains(newCortexKmer)) {
                                            kmerHeight = contigHeight/3;

                                            break;
                                        }
                                    }
                                }

                                for (byte[] e2 : SequenceUtils.generateSequencesWithEditDistance2(kmerBytes)) {
                                    CortexKmer newCortexKmer = new CortexKmer(e2);

                                    for (String altGene : geneKmers.keySet()) {
                                        if (geneKmers.get(altGene).contains(newCortexKmer)) {
                                            kmerHeight = contigHeight/2;

                                            break;
                                        }
                                    }
                                }

                                Color color = geneColors.get(geneName);

                                log.info("geneName={}, color={}", geneName, color);

                                stroke(color.getRGB(), 150.0f);
                                strokeCap(PROJECT);
                                line(xpos + k, ypos + 1, xpos + k, ypos + kmerHeight - 1);
                            }
                        }

                        fill(contigContainsKmerFromGene ? Color.RED.getRGB() : Color.BLACK.getRGB());
                        textSize(11);
                        textAlign(LEFT, TOP);
                        text(contig.hashCode(), marginX, ypos - 2);
                    }
                }
            }
        }
    }

    public void draw() {
        exit();
    }
}
