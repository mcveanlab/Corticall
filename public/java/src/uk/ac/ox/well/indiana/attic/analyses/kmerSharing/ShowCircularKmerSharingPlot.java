package uk.ac.ox.well.indiana.attic.analyses.kmerSharing;

import cern.colt.matrix.DoubleMatrix1D;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTreeMap;
import org.gicentre.utils.colour.ColourTable;
import uk.ac.ox.well.indiana.commands.Sketch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataFrame;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.statistics.misc.PCA;

import java.awt.*;
import java.io.PrintStream;
import java.util.*;

public class ShowCircularKmerSharingPlot extends Sketch {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="maxPrincipalCompoents", shortName="maxPC", doc="Maximum number of PCs to use when coloring kmer sharing plot")
    public Integer MAX_PRINCIPAL_COMPONENTS = 5;

    @Output
    public PrintStream out;

    public KmerSharingMatrix ksm = new KmerSharingMatrix();
    public IntervalTreeMap<TreeSet<String>> coordMapX = new IntervalTreeMap<TreeSet<String>>();
    public IntervalTreeMap<TreeSet<String>> coordMapY = new IntervalTreeMap<TreeSet<String>>();
    public HashMap<String, ScreenCoordinate> coords = new HashMap<String, ScreenCoordinate>();

    private class KmerSharingMatrix {
        public HashMap<String, HashMap<String, Integer>> matrix = new HashMap<String, HashMap<String, Integer>>();

        public void increment(String gene1, String gene2) {
            if (!matrix.containsKey(gene1)) { matrix.put(gene1, new HashMap<String, Integer>()); }
            if (!matrix.get(gene1).containsKey(gene2)) { matrix.get(gene1).put(gene2, 0); }

            matrix.get(gene1).put(gene2, matrix.get(gene1).get(gene2) + 1);
        }

        public int get(String gene1, String gene2) {
            return matrix.get(gene1).get(gene2);
        }

        public Set<String> keySet() {
            return matrix.keySet();
        }

        public Set<String> keySet(String gene1) {
            return matrix.get(gene1).keySet();
        }
    }

    private class ScreenCoordinate {
        public float x;
        public float y;

        public ScreenCoordinate(float x, float y) {
            this.x = x;
            this.y = y;
        }
    }

    private boolean isSharedKmer(CortexRecord cr) {
        int[] coverages = cr.getCoverages();

        boolean hasCoverage = false;
        int numColorsWithKmer = 0;
        boolean hasZeroOrUnitCoverageInColors = true;

        int totalCoverageInROI = 0;

        for (int color = 1; color < CORTEX_GRAPH.getNumColors(); color++) {
            hasCoverage |= (coverages[color] > 0);
            numColorsWithKmer += coverages[color] == 0 ? 0 : 1;
            hasZeroOrUnitCoverageInColors &= (coverages[color] <= 1);
            totalCoverageInROI += coverages[color];
        }

        boolean allCoverageIsInROI = (totalCoverageInROI == coverages[0]);

        return (hasCoverage && allCoverageIsInROI && numColorsWithKmer > 1 && hasZeroOrUnitCoverageInColors);
    }

    private HashMap<String, Integer> sampleColor = new HashMap<String, Integer>();
    private ColourTable ct;

    public void setup() {
        DataFrame<String, String, Float> d = new DataFrame<String, String, Float>(0.0f);

        int recordNum = 0;
        for (CortexRecord cr : CORTEX_GRAPH) {
            if (recordNum % (CORTEX_GRAPH.getNumRecords()/10) == 0) {
                log.info("processed {}/{} records", recordNum, CORTEX_GRAPH.getNumRecords());
            }
            recordNum++;

            if (isSharedKmer(cr)) {
                ArrayList<String> genesThatShareKmer = new ArrayList<String>();

                int[] coverages = cr.getCoverages();
                for (int color = 1; color < CORTEX_GRAPH.getNumColors(); color++) {
                    if (coverages[color] > 0) {
                        genesThatShareKmer.add(CORTEX_GRAPH.getColor(color).getSampleName());

                        d.set(cr.getKmerAsString(), CORTEX_GRAPH.getColor(color).getSampleName(), 1.0f);
                    }
                }

                for (int i = 0; i < genesThatShareKmer.size(); i++) {
                    for (int j = 0; j < genesThatShareKmer.size(); j++) {
                        if (i != j) {
                            String gene1 = genesThatShareKmer.get(i);
                            String gene2 = genesThatShareKmer.get(j);

                            ksm.increment(gene1, gene2);
                        }
                    }
                }
            }
        }

        log.info("performing PCA");
        PCA<String, String> pca = new PCA<String, String>(d);

        ct = ColourTable.getPresetColourTable(ColourTable.RD_YL_BU);

        log.info("minIndex={}, maxIndex={}", ct.getMinIndex(), ct.getMaxIndex());

        float colorIncrement = (Math.abs(ct.getMinIndex()) + Math.abs(ct.getMaxIndex())) / MAX_PRINCIPAL_COMPONENTS;

        TreeMap<Integer, TreeSet<String>> colorOrderMap = new TreeMap<Integer, TreeSet<String>>();

        for (int color = 1; color < CORTEX_GRAPH.getNumColors(); color++) {
            String sample = CORTEX_GRAPH.getColor(color).getSampleName();

            DoubleMatrix1D ev = pca.getEigenvector(sample);

            double maxPC = Math.abs(ev.get(0));
            int maxPCIndex = 0;
            for (int i = 1; i < MAX_PRINCIPAL_COMPONENTS; i++) {
                if (Math.abs(ev.get(i)) > maxPC) {
                    maxPC = Math.abs(ev.get(i));
                    maxPCIndex = i;
                }
            }

            float colorIndex = ct.getMinIndex() + maxPCIndex*colorIncrement;

            sampleColor.put(sample, ct.findColour(colorIndex));

            if (!colorOrderMap.containsKey(maxPCIndex)) {
                colorOrderMap.put(maxPCIndex, new TreeSet<String>());
            }
            colorOrderMap.get(maxPCIndex).add(sample);

            log.info("\tsample={}, color={}", sample, colorIndex);
        }

        //ArrayList<Integer> colorOrder = new ArrayList<Integer>();
        //for (Integer pcIndex : colorOrderMap.keySet()) {
            //TreeSet<String> samples
        //}

        size(displayHeight - 100, displayHeight - 100);

        drawBackground();
        drawLabels();
        drawLegend();
    }

    public void draw() {}

    private void drawLegend() {
        float colorIncrement = (Math.abs(ct.getMinIndex()) + Math.abs(ct.getMaxIndex())) / MAX_PRINCIPAL_COMPONENTS;

        int legendElementWidth = 30;

        for (int i = 0; i < MAX_PRINCIPAL_COMPONENTS; i++) {
            float colorIndex = ct.getMinIndex() + i*colorIncrement;

            fill(ct.findColour(colorIndex), 150.0f);
            noStroke();
            rect(-width/2 + 50 + legendElementWidth*i + legendElementWidth, height/2 - 100, legendElementWidth, legendElementWidth);

            fill(Color.BLACK.getRGB());
            textAlign(CENTER);
            text(String.format("PC%d", i + 1), -width / 2 + 50 + legendElementWidth * i + legendElementWidth + legendElementWidth / 2, height / 2 - 100 - legendElementWidth / 2);
        }
    }

    private void drawBackground() {
        translate(width/2, height/2);

        background(Color.WHITE.getRGB());

        strokeWeight(5);
        stroke(Color.BLACK.getRGB());
        fill(Color.WHITE.getRGB());

        float diameter = 0.8f*width;
        float radius = diameter/2;
        ellipse(0, 0, diameter, diameter);

        float angleDelta = radians(360.0f / (((float) CORTEX_GRAPH.getNumColors() - 1)));

        for (int color = 1; color < CORTEX_GRAPH.getNumColors(); color++) {
            float angle = ((float) color)*angleDelta;

            float xe = (radius * cos(angle));
            float ye = (radius * sin(angle));

            ellipse(xe, ye, 10, 10);

            String colorName = CORTEX_GRAPH.getColor(color).getSampleName();

            if (!coords.containsKey(colorName)) {
                coords.put(colorName, new ScreenCoordinate(xe, ye));

                Interval intervalX = new Interval("screen", (int) (xe - 5.0f), (int) (xe + 5.0f));
                Interval intervalY = new Interval("screen", (int) (ye - 5.0f), (int) (ye + 5.0f));

                if (!coordMapX.containsKey(intervalX)) {
                    coordMapX.put(intervalX, new TreeSet<String>());
                }
                coordMapX.get(intervalX).add(colorName);

                if (!coordMapY.containsKey(intervalY)) {
                    coordMapY.put(intervalY, new TreeSet<String>());
                }
                coordMapY.get(intervalY).add(colorName);
            }
        }
    }

    private void drawLabels() {
        float diameter = 0.8f*width;
        float radius = diameter/2;

        float angleDelta = radians(360.0f / (((float) CORTEX_GRAPH.getNumColors() - 1)));

        for (int color = 1; color < CORTEX_GRAPH.getNumColors(); color++) {
            rotate(angleDelta);

            fill(sampleColor.get(CORTEX_GRAPH.getColor(color).getSampleName()));
            textAlign(LEFT, CENTER);
            text(CORTEX_GRAPH.getColor(color).getSampleName(), radius + 10, 0);
        }
    }

    private Collection<String> intersectingElements(Collection<TreeSet<String>> s1, Collection<TreeSet<String>> s2) {
        TreeSet<String> genes1 = new TreeSet<String>();
        TreeSet<String> genes2 = new TreeSet<String>();

        for (TreeSet<String> g1 : s1) {
            genes1.addAll(g1);
        }

        for (TreeSet<String> g2 : s2) {
            genes2.addAll(g2);
        }

        genes1.retainAll(genes2);

        return genes1;
    }

    public void mousePressed() {
        drawBackground();

        int tx = mouseX - width/2;
        int ty = mouseY - height/2;

        Collection<TreeSet<String>> genesX = coordMapX.getOverlapping(new Interval("screen", tx, tx));
        Collection<TreeSet<String>> genesY = coordMapY.getOverlapping(new Interval("screen", ty, ty));

        Collection<String> genes = intersectingElements(genesX, genesY);

        strokeCap(SQUARE);
        noFill();

        for (String gene1 : genes) {
            for (String gene2 : ksm.keySet(gene1)) {
                int count = ksm.get(gene1, gene2);

                ScreenCoordinate sc1 = coords.get(gene1);
                ScreenCoordinate sc2 = coords.get(gene2);

                stroke(sampleColor.get(gene2), 150.0f);
                strokeWeight(1.2f*log(count) + 1);
                curve(sc1.x - 50, sc1.y - 100, sc1.x, sc1.y, sc2.x, sc2.y, sc2.x - 50, sc2.y - 100);
            }
        }

        drawLabels();
        drawLegend();
    }
}
