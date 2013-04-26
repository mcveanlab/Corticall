package uk.ac.ox.well.indiana.analyses.kmerSharing;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTreeMap;
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
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;

public class ShowKmerPlot extends Sketch {
    @Argument(fullName="reference", shortName="R", doc="Reference genome")
    public IndexedFastaSequenceFile FASTA;

    @Argument(fullName="gff", shortName="gff", doc="GFF of regions of interest")
    public GFF3 GFF;

    @Argument(fullName="genes", shortName="genes", doc="IDs of genes to evaluate")
    public ArrayList<String> GENES;

    @Argument(fullName="pca", shortName="pca", doc="PCA table", required=false)
    public File PCA;

    @Argument(fullName="pcaOn", shortName="pcaOn", doc="Is the PCA on genes or kmers?", required=false)
    public String PCA_ON = "KMERS";

    @Argument(fullName="kmerReferencePanel", shortName="krp", doc="Kmer reference panel to color", required=false)
    public File KMER_REFERENCE_PANEL;

    @Output
    public File out;

    private class GeneView {
        private GFF3Record gene;
        private ArrayList<GFF3Record> exons = new ArrayList<GFF3Record>();
        private IntervalTreeMap<Boolean> exonIntervals = new IntervalTreeMap<Boolean>();
        private String seq;
        private int kmerSize;
        private Color color;
        private ArrayList<GFF3Record> domains;
        private HashMap<String, Color> domainColors;

        public GeneView(Collection<GFF3Record> records, String seq, int kmerSize, Color color, ArrayList<GFF3Record> domains) {
            this.seq = seq;
            this.kmerSize = kmerSize;
            this.color = color;
            this.domains = domains;

            for (GFF3Record record : records) {
                if ("gene".equalsIgnoreCase(record.getType())) {
                    gene = record;
                } else if ("exon".equalsIgnoreCase(record.getType())) {
                    exons.add(record);

                    exonIntervals.put(new Interval(record.getSeqid(), record.getStart(), record.getEnd()), true);
                }
            }

            domainColors = new HashMap<String, Color>();

            domainColors.put("ATS", new Color(93, 46, 140));
            domainColors.put("CIDRa", new Color(5, 199, 242));
            domainColors.put("CIDRb", new Color(5, 199, 242));
            domainColors.put("CIDRd", new Color(5, 199, 242));
            domainColors.put("CIDRg", new Color(5, 199, 242));
            domainColors.put("CIDRpam", new Color(5, 199, 242));
            domainColors.put("DBLa", new Color(242, 29, 47));
            domainColors.put("DBLb", new Color(242, 29, 47));
            domainColors.put("DBLd", new Color(242, 29, 47));
            domainColors.put("DBLe", new Color(242, 29, 47));
            domainColors.put("DBLg", new Color(242, 29, 47));
            domainColors.put("DBLpam1", new Color(242, 29, 47));
            domainColors.put("DBLpam2", new Color(242, 29, 47));
            domainColors.put("DBLpam3", new Color(242, 29, 47));
            domainColors.put("DBLz", new Color(242, 29, 47));
            domainColors.put("NTS", new Color(242, 110, 34));
            domainColors.put("NTSpam", new Color(242, 110, 34));
        }

        public GFF3Record getGene() {
            return this.gene;
        }

        public boolean isExonic(int pos) {
            Interval interval;

            if (gene.getStrand() == GFF3Record.Strand.POSITIVE) {
                interval = new Interval(gene.getSeqid(), gene.getStart() + pos, gene.getStart() + pos);
            } else {
                interval = new Interval(gene.getSeqid(), gene.getEnd() - pos, gene.getEnd() - pos);
            }

            return exonIntervals.getOverlapping(interval).size() > 0;
        }

        public int getGeneLength() {
            return gene.getEnd() - gene.getStart();
        }

        public void display(int xpos0, int ypos, int exonHeight) {
            fill(0, 0, 0);
            textMode(SHAPE);
            textSize(33);
            textAlign(LEFT, CENTER);
            text(gene.getAttributes().get("ID"), 10, ypos);

            stroke(0, 0, 0);
            strokeWeight(1);
            line(xpos0, ypos, xpos0 + this.getGeneLength(), ypos);

            for (GFF3Record exon : exons) {
                int exonLength = exon.getEnd() - exon.getStart();
                int eypos0 = ypos - (exonHeight/2);

                int expos0;
                if (exon.getStrand() == GFF3Record.Strand.POSITIVE) {
                    expos0 = xpos0 + exon.getStart() - gene.getStart();
                } else {
                    expos0 = xpos0 + gene.getEnd() - exon.getEnd();
                }

                stroke(0, 0, 0);
                strokeWeight(1);
                fill(255, 255, 255);
                rect(expos0, eypos0, exonLength, exonHeight);
            }

            for (GFF3Record domain : domains) {
                int domainLength = domain.getEnd() - domain.getStart();

                int dxpos0 = xpos0 + domain.getStart() - gene.getStart();
                int dypos = ypos + (exonHeight/2) + 15;

                Color c = domainColors.get(domain.getAttributes().get("DOMAIN_TYPE"));

                stroke(c.getRGB());
                strokeWeight(8);
                line(dxpos0, dypos, dxpos0 + domainLength, dypos);
            }
        }
    }

    private class GeneViews extends ArrayList<GeneView> {
        private class KmerLoc {
            public int pos;
            public int row;
            public boolean isExonic;

            public KmerLoc(int pos, int row, boolean isExonic) {
                this.pos = pos;
                this.row = row;
                this.isExonic = isExonic;
            }
        }

        private class KmerMetaData {
            public Color color;
            public ArrayList<KmerLoc> kmerLocs = new ArrayList<KmerLoc>();
            public boolean display = false;

            public KmerMetaData(Color color) { this.color = color; }
            public boolean add(KmerLoc kmerLoc) { return kmerLocs.add(kmerLoc); }
        }

        private int horizontalMargin = 50;
        private int verticalMargin = 90;
        private int labelMargin = 260;
        private int exonHeight = 25;

        private int longestGeneLength = 0;
        private HashMap<String, KmerMetaData> kmers = new HashMap<String, KmerMetaData>();

        @Override
        public boolean add(GeneView e) {
            int geneLength = e.getGeneLength();

            if (geneLength > longestGeneLength) {
                longestGeneLength = geneLength;
            }

            for (int i = 0; i < e.seq.length() - e.kmerSize; i++) {
                String kmer = SequenceUtils.alphanumericallyLowestOrientation(e.seq.substring(i, i + e.kmerSize));

                if (!kmers.containsKey(kmer)) {
                    kmers.put(kmer, new KmerMetaData(e.color));
                }

                if (e.getGene().getStrand() == GFF3Record.Strand.POSITIVE) {
                    kmers.get(kmer).add(new KmerLoc(i, this.size(), e.isExonic(i)));
                } else {
                    kmers.get(kmer).add(new KmerLoc(e.seq.length() - e.kmerSize - i, this.size(), e.isExonic(e.seq.length() - e.kmerSize - i)));
                }
            }

            return super.add(e);
        }

        public int getDisplayWidth() {
            return labelMargin + horizontalMargin + longestGeneLength + horizontalMargin;
        }

        public int getDisplayHeight() {
            return this.size() * verticalMargin;
        }

        public void display() {
            int xpos0 = labelMargin + (horizontalMargin/2);

            for (int geneIndex = 0; geneIndex < this.size(); geneIndex++) {
                int ypos  = (geneIndex*verticalMargin) + (verticalMargin/2);

                this.get(geneIndex).display(xpos0, ypos, exonHeight);
            }

            strokeCap(PROJECT);
            strokeWeight(1);

            for (String kmer : kmers.keySet()) {
                KmerMetaData kmd = kmers.get(kmer);

                if (kmd.display) {
                    stroke(kmd.color.getRed(), kmd.color.getGreen(), kmd.color.getBlue());

                    for (KmerLoc kmerLoc : kmd.kmerLocs) {
                        int kxpos = kmerLoc.pos;
                        int kypos = (kmerLoc.row*verticalMargin) + (verticalMargin/2);

                        if (kmerLoc.isExonic) {
                            line(xpos0 + kxpos, kypos - (exonHeight/2) + 1, xpos0 + kxpos, kypos + (exonHeight/2) - 1);
                        } else {
                            line(xpos0 + kxpos, kypos - (exonHeight/5) + 1, xpos0 + kxpos, kypos + (exonHeight/5) - 1);
                        }
                    }
                }
            }
        }
    }

    private class Legend {
        private int xpos = 0;
        private int ypos = 0;

        private int verticalMargin = 40;
        private int lineLength = 20;
        private int lineAndLabelMargin = 10;

        ArrayList<String> labels = new ArrayList<String>();
        ArrayList<Color> colors = new ArrayList<Color>();

        public Legend(int xpos, int ypos) {
            this.xpos = xpos;
            this.ypos = ypos;
        }

        public void addElement(String label, Color color) {
            labels.add(label);
            colors.add(color);
        }

        public void display() {
            for (int i = 0; i < labels.size(); i++) {
                int ypos1 = ypos + i*verticalMargin;

                strokeWeight(7);
                strokeCap(PROJECT);
                stroke(colors.get(i).getRed(), colors.get(i).getGreen(), colors.get(i).getBlue());
                line(xpos, ypos1, xpos + lineLength, ypos1);

                fill(0, 0, 0);
                textMode(SHAPE);
                textSize(32);
                textAlign(LEFT, CENTER);
                text(labels.get(i), xpos + lineLength + lineAndLabelMargin, ypos1);
            }
        }

        public int getDisplayWidth() {
            return (lineLength + lineAndLabelMargin + 500);
        }

        public int getDisplayHeight() {
            return labels.size()*verticalMargin;
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

    private HashSet<String> loadKmerReferencePanel() {
        HashSet<String> krp = new HashSet<String>();

        TableReader tr = new TableReader(KMER_REFERENCE_PANEL);

        for (HashMap<String, String> entry : tr) {
            krp.add(entry.get("kmer"));
        }

        return krp;
    }

    private int getKmerSize(HashSet<String> krp) {
        return krp.iterator().next().length();
    }

    private GeneViews geneViews = null;
    private Legend legend = null;
    private Legend domainLegend = null;

    public void setup() {
        HashSet<String> krp = loadKmerReferencePanel();
        int kmerSize = getKmerSize(krp);

        HashMap<String, HashMap<String, String>> pcaTable = new HashMap<String, HashMap<String, String>>();
        if (PCA != null) {
            TableReader tr = new TableReader(PCA);
            for (HashMap<String, String> fields : tr) {
                pcaTable.put(fields.get(""), fields);
            }
        }

        Color[] colors = generateColors(GENES.size());

        String[] columns = new String[] { "PC1", "PC2", "PC3", "PC4", "PC5", "PC6" };

        if (PCA_ON.equalsIgnoreCase("genes")) {
            Color[] newcolors = generateColors(columns.length);

            for (int geneIndex = 0; geneIndex < GENES.size(); geneIndex++) {
                String gene = GENES.get(geneIndex);

                float maxValue = 0.0f;
                int maxIndex = 0;
                for (int i = 0; i < columns.length; i++) {
                    float value = Math.abs(Float.parseFloat(pcaTable.get(gene).get(columns[i])));

                    if (value > maxValue) {
                        maxValue = value;
                        maxIndex = i;
                    }
                }

                colors[geneIndex] = newcolors[maxIndex];
            }
        }

        geneViews = new GeneViews();

        int geneIndex = 0;
        for (String gene : GENES) {
            GFF3Record record = GFF.getRecord(gene);

            String seq = new String(FASTA.getSubsequenceAt(record.getSeqid(), record.getStart(), record.getEnd()).getBases());

            Collection<GFF3Record> records = GFF.getOverlapping(record);

            ArrayList<GFF3Record> domains = new ArrayList<GFF3Record>();

            for (GFF3Record r : records) {
                if ("domain".equalsIgnoreCase(r.getType())) {
                    domains.add(r);
                }
            }

            geneViews.add(new GeneView(GFF.getContained(record), seq, kmerSize, colors[geneIndex], domains));

            geneIndex++;
        }

        legend = new Legend(geneViews.getDisplayWidth(), geneViews.verticalMargin);
        if (PCA != null) {
            for (int i = 0; i < columns.length; i++) {
                legend.addElement(columns[i], colors[i]);
            }
        } else {
            for (int i = 0; i < GENES.size(); i++) {
                legend.addElement(GENES.get(i), colors[i]);
            }
        }

        domainLegend = new Legend(geneViews.getDisplayWidth(), legend.getDisplayHeight() + 5*geneViews.verticalMargin);

        domainLegend.addElement("ATS", new Color(93, 46, 140));
        domainLegend.addElement("CIDR", new Color(5, 199, 242));
        domainLegend.addElement("DBL", new Color(242, 29, 47));
        domainLegend.addElement("NTS", new Color(242, 110, 34));

        size(geneViews.getDisplayWidth() + legend.getDisplayWidth(), geneViews.getDisplayHeight(), PGraphicsPDF.PDF, out.getAbsolutePath());
        noLoop();

        if (PCA != null) {
            colors = generateColors(10);
        }

        int krpIndex = 0;
        for (String kmer : krp) {
            if (geneViews.kmers.containsKey(kmer)) {
                geneViews.kmers.get(kmer).display = true;

                if (PCA != null && PCA_ON.equalsIgnoreCase("kmers") && pcaTable.containsKey(kmer)) {
                    float maxValue = 0.0f;
                    int maxIndex = 0;
                    for (int i = 0; i < columns.length; i++) {
                        float value = Math.abs(Float.parseFloat(pcaTable.get(kmer).get(columns[i])));

                        if (value > maxValue) {
                            maxValue = value;
                            maxIndex = i;
                        }
                    }

                    geneViews.kmers.get(kmer).color = colors[maxIndex];
                }
            }

            if (krpIndex % (krp.size() / 10) == 0) {
                log.info("Processed {}/{} records", krpIndex, krp.size());
            }
            krpIndex++;
        }

        log.info("Gene display: {}x{}", geneViews.getDisplayWidth(), geneViews.getDisplayHeight());
        log.info("Legend display: {}x{}", legend.getDisplayWidth(), legend.getDisplayHeight());
    }

    public void draw() {
        geneViews.display();
        legend.display();
        domainLegend.display();

        exit();
    }
}
