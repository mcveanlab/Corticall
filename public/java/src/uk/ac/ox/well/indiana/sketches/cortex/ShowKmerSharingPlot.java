package uk.ac.ox.well.indiana.sketches.cortex;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTree;
import net.sf.picard.util.IntervalTreeMap;
import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlContext;
import org.apache.commons.jexl2.MapContext;
import org.geotools.brewer.color.ColorBrewer;
import org.geotools.brewer.color.StyleGenerator;
import uk.ac.ox.well.indiana.sketches.Sketch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.gff.*;

import java.awt.*;
import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

import processing.pdf.PGraphicsPDF;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

public class ShowKmerSharingPlot extends Sketch {
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

    private class Gene {
        public IntervalTreeMap<Integer> exonIntervals = new IntervalTreeMap<Integer>();

        public Gene(Collection<GFF3Record> records) {
            for (GFF3Record record : records) {
                if ("exon".equalsIgnoreCase(record.getType())) {
                    Interval interval = new Interval(record.getSeqid(), record.getStart(), record.getEnd());

                    exonIntervals.put(interval, null);
                }
            }
        }

        public boolean isExonic(String seq, int start) {
            Interval interval = new Interval(seq, start, start);

            return (exonIntervals.getOverlapping(interval).size() > 0);
        }
    }

    private class KmerMetaData {
        public Color color = null;
        public ArrayList<Integer> xpos = new ArrayList<Integer>();
        public ArrayList<Integer> ypos = new ArrayList<Integer>();
        public ArrayList<Boolean> isExonic = new ArrayList<Boolean>();

        public void addEntry(Color color, int xpos, int ypos, boolean isExonic) {
            if (this.color == null) {
                this.color = color;
            }

            this.xpos.add(xpos);
            this.ypos.add(ypos);
            this.isExonic.add(isExonic);
        }
    }

    // Taken from http://stackoverflow.com/questions/223971/how-to-generate-spectrum-color-palettes/223981#223981
    public Color[] generateColors(int n) {
        Color[] cols = new Color[n];

        for(int i = 0; i < n; i++) {
            cols[i] = Color.getHSBColor((float) i / (float) n, 0.85f, 1.0f);
        }
        return cols;
    }

    public void setup() {
        HashMap<String, Collection<GFF3Record>> allGeneRecords = new HashMap<String, Collection<GFF3Record>>();

        HashMap<String, Gene> genes = new HashMap<String, Gene>();

        // Figure out longest gene
        int longestGeneLength = 0;

        for (String gene : GENES) {
            GFF3Record record = GFF.getRecord(gene);
            Collection<GFF3Record> records = GFF.getContained(record);

            genes.put(gene, new Gene(records));

            allGeneRecords.put(gene, records);

            int geneLength = record.getEnd() - record.getStart();
            if (geneLength > longestGeneLength) {
                longestGeneLength = geneLength;
            }
        }

        // Set display properties
        int verticalMargin = 70;
        int horizontalMargin = 50;
        int exonHeight = 50;

        int windowWidth = longestGeneLength + horizontalMargin;
        int windowHeight = GENES.size() * verticalMargin;

        size(windowWidth, windowHeight, PGraphicsPDF.PDF, out.getAbsolutePath());

        // Load up kmers
        HashMap<String, KmerMetaData> kmerData = new HashMap<String, KmerMetaData>();

        Color[] colors = generateColors(allGeneRecords.size());

        int geneIndex = 0;
        for (String gene : allGeneRecords.keySet()) {
            Gene theGene = genes.get(gene);

            GFF3Record geneRecord = null;

            for (GFF3Record record : allGeneRecords.get(gene)) {
                if ("gene".equalsIgnoreCase(record.getType())) {
                    geneRecord = record;

                    int geneLength = geneRecord.getEnd() - geneRecord.getStart();

                    int xpos0 = (horizontalMargin/2);
                    int xpos1 = (horizontalMargin/2) + geneLength;
                    int ypos = (geneIndex*verticalMargin)+(verticalMargin/2);

                    line(xpos0, ypos, xpos1, ypos);

                    String seq = new String(FASTA.getSubsequenceAt(record.getSeqid(), record.getStart(), record.getEnd()).getBases());
                    for (int i = 0; i < seq.length() - CORTEX_GRAPH.getKmerSize(); i++) {
                        String kmer = SequenceUtils.getAlphanumericallyLowestOrientation(seq.substring(i, i + CORTEX_GRAPH.getKmerSize()));

                        if (!kmerData.containsKey(kmer)) {
                            kmerData.put(kmer, new KmerMetaData());
                        }

                        kmerData.get(kmer).addEntry(colors[geneIndex], xpos0 + i, ypos, theGene.isExonic(record.getSeqid(), record.getStart() + i));
                    }
                }
            }

            if (geneRecord != null) {
                for (GFF3Record record : allGeneRecords.get(gene)) {
                    if ("exon".equalsIgnoreCase(record.getType())) {
                        int exonLength = record.getEnd() - record.getStart();

                        int xpos0 = (horizontalMargin/2) + (record.getStart() - geneRecord.getStart());
                        int ypos0 = (geneIndex*verticalMargin)+(verticalMargin/2)-(exonHeight/2);
                        int width = exonLength;
                        int height = exonHeight;

                        rect(xpos0, ypos0, width, height);
                    }
                }
            }

            geneIndex++;
        }

        int recordNum = 0;
        for (CortexRecord cr : CORTEX_GRAPH) {
            int[] coverages = cr.getCoverages();

            boolean hasCoverage = false;
            int numColorsWithKmer = 0;
            boolean hasZeroOrUnitCoverageInColors = true;

            for (int color : COLORS) {
                hasCoverage |= (coverages[color] > 0);
                numColorsWithKmer += coverages[color] == 0 ? 0 : 1;
                hasZeroOrUnitCoverageInColors &= (coverages[color] <= 1);
            }

            if (hasCoverage && numColorsWithKmer > 1 && hasZeroOrUnitCoverageInColors && kmerData.containsKey(cr.getKmerString())) {
                KmerMetaData kmd = kmerData.get(cr.getKmerString());

                for (int i = 0; i < kmd.xpos.size(); i++) {
                    int xpos = kmd.xpos.get(i);
                    int ypos = kmd.ypos.get(i);
                    boolean isExonic = kmd.isExonic.get(i);

                    stroke(kmd.color.getRed(), kmd.color.getGreen(), kmd.color.getBlue());

                    if (isExonic) {
                        line(xpos, ypos-(exonHeight/2)+1, xpos, ypos+(exonHeight/2)-1);
                    } else {
                        line(xpos, ypos-(exonHeight/4)+1, xpos, ypos+(exonHeight/4)-1);
                    }
                }
            }

            recordNum++;

            if (recordNum % (CORTEX_GRAPH.getNumRecords()/10) == 0) {
                log.info("processed {}/{} records", recordNum, CORTEX_GRAPH.getNumRecords());
            }
        }

        exit();
    }
}
