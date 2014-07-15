package uk.ac.ox.well.indiana.attic.analyses.nahr;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import processing.pdf.PGraphicsPDF;
import uk.ac.ox.well.indiana.commands.Sketch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;

import java.awt.*;
import java.io.File;
import java.util.*;

public class PaintParentalContributionsOnContigs extends Sketch {
    @Argument(fullName="contigs", shortName="c", doc="Aligned contigs (BAM)")
    public SAMFileReader CONTIGS;

    @Argument(fullName="parentGraph", shortName="pg", doc="Parent graph")
    public TreeMap<String, CortexGraph> PARENT_GRAPH;

    @Argument(fullName="parentColor", shortName="pc", doc="Parent color")
    public TreeMap<String, Color> PARENT_COLOR;

    @Output
    public File out;

    private final int marginXRight = 300;
    private final int marginLabel = 285;
    private final int marginX = 10;
    private final int marginY = 10;
    private final int marginTitle = 30;
    private final int marginContig = 5;
    private final int contigHeight = 10;

    private int numContigs = 0;
    private int maxContigLength = 0;
    private int kmerSize = 0;

    private final Color sharedColor = Color.WHITE;
    private final Color absenceColor = Color.LIGHT_GRAY;

    private Map<CortexKmer, Color> contigKmerColors = new HashMap<CortexKmer, Color>();
    private Set<SAMRecord> contigs = new TreeSet<SAMRecord>(new Comparator<SAMRecord>() {
        @Override
        public int compare(SAMRecord o1, SAMRecord o2) {
            if (o1.getReadLength() == o2.getReadLength()) { return 0; }

            return o1.getReadLength() > o2.getReadLength() ? 1 : -1;
        }
    });

    public void initialize() {
        kmerSize = PARENT_GRAPH.values().iterator().next().getKmerSize();

        log.info("Loading contigs...");
        for (SAMRecord contig : CONTIGS) {
            numContigs++;

            if (contig.getReadLength() > maxContigLength) {
                maxContigLength = contig.getReadLength();
            }

            String seq = contig.getReadString();
            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                CortexKmer kmer = new CortexKmer(seq.substring(i, i + kmerSize));

                contigKmerColors.put(kmer, null);
            }

            contigs.add(contig);
        }

        log.info("Assigning colors from parental graphs...");
        for (String parentName : PARENT_GRAPH.keySet()) {
            log.info("  assigning kmers from {} the color {}", parentName, PARENT_COLOR.get(parentName));

            CortexGraph parentGraph = PARENT_GRAPH.get(parentName);

            int recordsSeen = 0;
            for (CortexRecord cr : parentGraph) {
                CortexKmer ck = cr.getKmer();

                if (contigKmerColors.containsKey(ck)) {
                    if (contigKmerColors.get(ck) == null) {
                        contigKmerColors.put(ck, PARENT_COLOR.get(parentName));
                    } else {
                        contigKmerColors.put(ck, sharedColor);
                    }
                }

                if (recordsSeen % (parentGraph.getNumRecords() / 5) == 0) {
                    log.info("    processed {}/{} (~{}%) records", recordsSeen, parentGraph.getNumRecords(), 100*recordsSeen/parentGraph.getNumRecords());
                }
                recordsSeen++;
            }
        }
    }

    public void setup() {
        size(marginLabel + 2*marginX + maxContigLength + marginXRight, marginTitle + 2*marginY + (marginContig + contigHeight)*(numContigs + 10), PGraphicsPDF.PDF, out.getAbsolutePath());

        background(Color.WHITE.getRGB());
        strokeCap(PROJECT);

        int index = 0;
        for (SAMRecord contig : contigs) {
            int xpos = marginLabel + marginX;
            int ypos = marginTitle + marginY + index*(marginContig + contigHeight);

            String label = (index+1) + ". " + contig.getReadGroup().getSample() + " " + contig.getReadGroup().getId() + " " + contig.getReadName();

            textSize(11);
            textAlign(LEFT, TOP);
            fill(Color.BLACK.getRGB());
            text(label, marginX + 30, ypos - 2);

            String contigSeq = contig.getReadString();

            Map<String, Integer> kmerCounts = new LinkedHashMap<String, Integer>();
            for (String sampleName : PARENT_COLOR.keySet()) {
                kmerCounts.put(sampleName, 0);
            }
            kmerCounts.put("shared", 0);
            kmerCounts.put("absent", 0);

            Color[] colors = new Color[contigSeq.length()];

            for (int i = 0; i <= contigSeq.length() - kmerSize; i++) {
                CortexKmer kmer = new CortexKmer(contigSeq.substring(i, i + kmerSize));

                if (contigKmerColors.containsKey(kmer) && contigKmerColors.get(kmer) != null && !contigKmerColors.get(kmer).equals(sharedColor)) {
                    Color color = contigKmerColors.get(kmer);

                    for (int j = i; j < i + kmerSize; j++) {
                        if (colors[j] == null) {
                            colors[j] = color;
                        }
                    }

                    for (String sampleName : PARENT_COLOR.keySet()) {
                        if (contigKmerColors.get(kmer).equals(PARENT_COLOR.get(sampleName))) {
                            kmerCounts.put(sampleName, kmerCounts.get(sampleName) + 1);
                        }
                    }
                }
            }

            for (int i = 0; i <= contigSeq.length() - kmerSize; i++) {
                CortexKmer kmer = new CortexKmer(contigSeq.substring(i, i + kmerSize));

                if (contigKmerColors.containsKey(kmer)) {
                    if (contigKmerColors.get(kmer) == null) {
                        for (int j = i; j < i + kmerSize; j++) {
                            if (colors[j] == null) {
                                colors[j] = absenceColor;
                            }
                        }

                        kmerCounts.put("absent", kmerCounts.get("absent") + 1);
                    } else if (contigKmerColors.get(kmer).equals(sharedColor)) {
                        // draw nothing

                        kmerCounts.put("shared", kmerCounts.get("shared") + 1);
                    }
                }
            }

            for (int i = 0; i < contigSeq.length(); i++) {
                if (colors[i] != null) {
                    stroke(colors[i].getRGB(), 255.0f);
                    line(xpos + i, ypos, xpos + i, ypos + contigHeight);
                }
            }

            stroke(Color.BLACK.getRGB());
            noFill();
            rect(xpos, ypos, contig.getReadLength(), contigHeight);

            String contigLabel = Joiner.on("; ").withKeyValueSeparator("=").join(kmerCounts);
            textSize(11);
            textAlign(LEFT, TOP);
            fill(Color.BLACK.getRGB());
            text(contigLabel, xpos + contig.getReadLength() + 5, ypos - 2);

            index++;
        }
    }

    public void draw() {
        exit();
    }
}
