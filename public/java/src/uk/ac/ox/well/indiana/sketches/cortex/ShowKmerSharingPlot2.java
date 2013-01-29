package uk.ac.ox.well.indiana.sketches.cortex;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTreeMap;
import processing.pdf.PGraphicsPDF;
import uk.ac.ox.well.indiana.sketches.Sketch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.awt.*;
import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

public class ShowKmerSharingPlot2 extends Sketch {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="colors", shortName="c", doc="Colors to process")
    public ArrayList<Integer> COLORS;

    @Argument(fullName="reference", shortName="R", doc="Reference genome")
    public IndexedFastaSequenceFile FASTA;

    @Argument(fullName="gff", shortName="gff", doc="GFF of regions of interest")
    public GFF3 GFF;

    @Argument(fullName="genes", shortName="genes", doc="IDs of genes to evaluate")
    public ArrayList<String> GENES;

    @Output
    public File out;

    private class GeneView {
        private GFF3Record gene;
        private ArrayList<GFF3Record> exons = new ArrayList<GFF3Record>();
        private IntervalTreeMap<Boolean> exonIntervals = new IntervalTreeMap<Boolean>();
        private String seq;
        private int kmerSize;
        private Color color;

        public GeneView(Collection<GFF3Record> records, String seq, int kmerSize, Color color) {
            this.seq = seq;
            this.kmerSize = kmerSize;
            this.color = color;

            for (GFF3Record record : records) {
                if ("gene".equalsIgnoreCase(record.getType())) {
                    gene = record;
                } else if ("exon".equalsIgnoreCase(record.getType())) {
                    exons.add(record);

                    exonIntervals.put(new Interval(record.getSeqid(), record.getStart(), record.getEnd()), true);
                }
            }
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
            textSize(32);
            text(gene.getAttributes().get("ID"), 10, ypos);

            line(xpos0, ypos, xpos0 + this.getGeneLength(), ypos);

            fill(255, 255, 255);

            for (GFF3Record exon : exons) {
                int exonLength = exon.getEnd() - exon.getStart();
                int eypos0 = ypos - (exonHeight/2);

                int expos0;
                if (exon.getStrand() == GFF3Record.Strand.POSITIVE) {
                    expos0 = xpos0 + exon.getStart() - gene.getStart();
                } else {
                    expos0 = xpos0 + gene.getEnd() - exon.getEnd();
                }

                rect(expos0, eypos0, exonLength, exonHeight);
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
        private int verticalMargin = 70;
        private int labelMargin = 260;
        private int exonHeight = 50;

        private int longestGeneLength = 0;
        private HashMap<String, KmerMetaData> kmers = new HashMap<String, KmerMetaData>();

        @Override
        public boolean add(GeneView e) {
            int geneLength = e.getGeneLength();

            if (geneLength > longestGeneLength) {
                longestGeneLength = geneLength;
            }

            for (int i = 0; i < e.seq.length() - e.kmerSize; i++) {
                String kmer = SequenceUtils.getAlphanumericallyLowestOrientation(e.seq.substring(i, i + e.kmerSize));

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

    // Taken from http://stackoverflow.com/questions/223971/how-to-generate-spectrum-color-palettes/223981#223981
    private Color[] generateColors(int n) {
        Color[] cols = new Color[n];

        for(int i = 0; i < n; i++) {
            cols[i] = Color.getHSBColor((float) i / (float) n, 0.85f, 1.0f);
        }

        return cols;
    }

    private GeneViews geneViews = new GeneViews();

    public void setup() {
        Color[] colors = generateColors(GENES.size());

        int geneIndex = 0;
        for (String gene : GENES) {
            GFF3Record record = GFF.getRecord(gene);

            String seq = new String(FASTA.getSubsequenceAt(record.getSeqid(), record.getStart(), record.getEnd()).getBases());

            geneViews.add(new GeneView(GFF.getContained(record), seq, CORTEX_GRAPH.getKmerSize(), colors[geneIndex]));

            geneIndex++;
        }

        int crindex = 0;
        for (CortexRecord cr : CORTEX_GRAPH) {
            int[] coverages = cr.getCoverages();

            boolean hasCoverage = false;
            int numColorsWithKmer = 0;
            boolean hasZeroOrUnitCoverageInColors = true;

            int totalCoverageInROI = 0;

            for (int color : COLORS) {
                hasCoverage |= (coverages[color] > 0);
                numColorsWithKmer += coverages[color] == 0 ? 0 : 1;
                hasZeroOrUnitCoverageInColors &= (coverages[color] <= 1);
                totalCoverageInROI += coverages[color];
            }

            boolean allCoverageIsInROI = (totalCoverageInROI == coverages[0]);

            if (hasCoverage && allCoverageIsInROI && numColorsWithKmer > 1 && hasZeroOrUnitCoverageInColors && geneViews.kmers.containsKey(cr.getKmerString())) {
                geneViews.kmers.get(cr.getKmerString()).display = true;
            }

            if (crindex % (CORTEX_GRAPH.getNumRecords() / 10) == 0) {
                log.info("Processed {}/{} records", crindex, CORTEX_GRAPH.getNumRecords());
            }
            crindex++;
        }

        size(geneViews.getDisplayWidth(), geneViews.getDisplayHeight(), PGraphicsPDF.PDF, out.getAbsolutePath());

        geneViews.display();

        exit();
    }
}
