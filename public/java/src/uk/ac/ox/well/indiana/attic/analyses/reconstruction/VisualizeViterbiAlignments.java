package uk.ac.ox.well.indiana.attic.analyses.reconstruction;

import com.google.common.base.Joiner;
import processing.core.PFont;
import processing.pdf.PGraphicsPDF;
import uk.ac.ox.well.indiana.commands.Sketch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class VisualizeViterbiAlignments extends Sketch {
    @Argument(fullName="alignment", shortName="a", doc="Alignment file")
    public File ALIGNMENT;

    @Argument(fullName="kmerReferencePanel", shortName="krp", doc="Kmer reference panel")
    public File KMER_REFERENCE_PANEL;

    @Output
    public File out;

    private class AlignmentEntry {
        public List<String> lines = new ArrayList<String>();

        public String toString() {
            return Joiner.on("\n").join(lines) + "\n";
        }
    }

    private List<AlignmentEntry> aes = new ArrayList<AlignmentEntry>();

    private final int marginTitle = 30;
    private final int marginX = 10;
    private final int marginY = 10;
    private final int contigHeight = 10;

    private int baseWidth = 5;
    private int longestLine = 0;
    private int maxLines = 0;

    private int sketchWidth = 0;
    private int sketchHeight = 0;
    private PFont pfont;

    private void loadAlignments() {
        try {
            AlignmentEntry ae = null;

            BufferedReader br = new BufferedReader(new FileReader(ALIGNMENT));
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#") || line.isEmpty()) {
                    // do nothing
                } else if (line.startsWith("Target")) {
                    if (ae != null) {
                        aes.add(ae);
                    }

                    ae = new AlignmentEntry();
                    ae.lines.add(line.replaceAll("\t", "    "));
                } else if (line.startsWith("Parameters")) {
                    aes.add(ae);

                    break;
                } else {
                    if (ae != null) {
                        ae.lines.add(line.replaceAll("\t", "    "));
                    }
                }
            }
        } catch (FileNotFoundException e) {
            throw new RuntimeException("File '" + ALIGNMENT.getAbsolutePath() + "' not found: " + e);
        } catch (IOException e) {
            throw new RuntimeException("Error reading '" + ALIGNMENT.getAbsolutePath() + "' not found: " + e);
        }
    }

    public void initialize() {
        loadAlignments();

        for (AlignmentEntry ae : aes) {
            for (String line : ae.lines) {
                int length = line.length();

                if (length > longestLine) {
                    longestLine = length;
                }
            }

            if (ae.lines.size() > maxLines) {
                maxLines = ae.lines.size();
            }
        }

        pfont = createFont("Monaco", 10);
        baseWidth = (int) textWidth('W');

        sketchWidth = marginX + (longestLine*baseWidth) + marginX;
        sketchHeight = marginY + marginTitle + (contigHeight*maxLines) + marginY;

        log.info("longestLine: {}", longestLine);
        log.info("maxLines: {}", maxLines);
        log.info("sketch dimensions: {}x{}", sketchWidth, sketchHeight);
    }

    public void setup() {
        size(sketchWidth, sketchHeight, PGraphicsPDF.PDF, out.getAbsolutePath());

        PGraphicsPDF pdf = (PGraphicsPDF) g;

        for (int i = 0; i < aes.size(); i++) {
            if (i > 0) { pdf.nextPage(); }

            background(204, 204, 204);

            AlignmentEntry ae = aes.get(i);

            for (int j = 0; j < ae.lines.size(); j++) {
                fill(Color.BLACK.getRGB());
                textFont(pfont);
                text(ae.lines.get(j), marginX, marginY + (j*contigHeight));
            }
        }

        exit();
    }
}
