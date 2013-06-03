package uk.ac.ox.well.indiana.analyses.reconstruction;

import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlContext;
import org.apache.commons.jexl2.MapContext;
import processing.pdf.PGraphicsPDF;
import uk.ac.ox.well.indiana.sketches.Sketch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

import java.awt.*;
import java.io.File;
import java.io.PrintStream;
import java.util.*;
import java.util.List;

public class VisualizeMultiPanelContigs extends Sketch {
    @Argument(fullName="contigTable", shortName="ct", doc="Contig table")
    public File CONTIG_TABLE;

    @Argument(fullName="onlyTheseGenes", shortName="g", doc="Process only these genes", required=false)
    public HashSet<String> ONLY_THESE_GENES;

    @Argument(fullName="onlyTheseSamples", shortName="s", doc="Process only these samples", required=false)
    public HashSet<String> ONLY_THESE_SAMPLES;

    @Output
    public File out;

    //          sample      gene    entry
    private Map<String, Map<String, List<Map<String, String>>>> data;

    //          gene    color
    private Map<String, Color> geneColors;

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
            String contig = te.get("contig");

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

        for (String sample : data.keySet()) {
            for (String gene : data.get(sample).keySet()) {
                if ( (ONLY_THESE_SAMPLES == null || ONLY_THESE_SAMPLES.contains(sample)) && (ONLY_THESE_GENES == null || ONLY_THESE_GENES.contains(gene)) ) {
                    if (data.get(sample).get(gene).size() > maxContigs) {
                        maxContigs = data.get(sample).get(gene).size();
                    }

                    for (Map<String, String> te : data.get(sample).get(gene)) {
                        String contig = te.get("contig");

                        if (contig.length() > longestContig) {
                            longestContig = contig.length();
                        }
                    }
                }
            }
        }

        Color[] colors = generateColors(geneColors.size());

        int colorIndex = 0;
        for (String gene : geneColors.keySet()) {
            geneColors.put(gene, colors[colorIndex]);

            colorIndex++;
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
        size(marginLabel + 2*marginX + longestContig, marginTitle + 2*marginY + (marginContig + contigHeight)*maxContigs, PGraphicsPDF.PDF, out.getAbsolutePath());

        PGraphicsPDF pdf = (PGraphicsPDF) g;

        log.info("{}x{}", getWidth(), getHeight());

        int index = 0;
        for (String sample : data.keySet()) {
            for (String gene : data.get(sample).keySet()) {
                if ( (ONLY_THESE_SAMPLES == null || ONLY_THESE_SAMPLES.contains(sample)) && (ONLY_THESE_GENES == null || ONLY_THESE_GENES.contains(gene)) ) {
                    if (index > 0) { pdf.nextPage(); }
                    index++;

                    log.info("{} {} {}", sample, gene, data.get(sample).get(gene).size());

                    background(204, 204, 204);

                    stroke(geneColors.get(gene).getRGB());
                    fill(geneColors.get(gene).getRGB());
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
                        textSize(11);
                        textAlign(LEFT, TOP);
                        text(contig.hashCode(), marginX, ypos - 2);

                        for (int k = 0; k <= contig.length() - kmerSize; k++) {
                            CortexKmer kmer = new CortexKmer(contig.substring(k, k+kmerSize));

                            if (kmersAndGenes.containsKey(kmer)) {
                                String geneName = kmersAndGenes.get(kmer);
                                Color color = geneColors.get(geneName);

                                stroke(color.getRGB());
                                strokeCap(PROJECT);
                                line(xpos + k, ypos + 1, xpos + k, ypos + contigHeight - 1);
                            }
                        }
                    }
                }
            }
        }
    }

    public void draw() {
        exit();
    }
}
