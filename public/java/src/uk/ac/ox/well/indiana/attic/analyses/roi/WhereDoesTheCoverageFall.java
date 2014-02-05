package uk.ac.ox.well.indiana.attic.analyses.roi;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import processing.core.PFont;
import processing.pdf.PGraphicsPDF;
import uk.ac.ox.well.indiana.commands.Sketch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;

import java.awt.*;
import java.io.File;
import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.text.DecimalFormat;
import java.util.*;

/**
 * This sketch shows where a specified list of kmers appears on the given genome (and if provided, the genes in which they fall)
 */
public class WhereDoesTheCoverageFall extends Sketch {
    @Argument(fullName="reference", shortName="R", doc="Reference sequence")
    public FastaSequenceFile REFERENCE;

    @Argument(fullName="gff", shortName="gff", doc="Gene annotations", required=false)
    public GFF3 GFF;

    @Argument(fullName="bam", shortName="B", doc="BAM file")
    public SAMFileReader BAM;

    @Argument(fullName="maxCoverage", shortName="mc", doc="Maximum coverage")
    public Integer MAX_COVERAGE = 300;

    @Argument(fullName="WINDOW", shortName="w", doc="Window size")
    public Integer WINDOW = 1000;

    @Output
    public File out;

    private class Locus {
        public String contig;
        public int pos;

        public Locus(String contig, int pos) {
            this.contig = contig;
            this.pos = pos;
        }
    }

    private Map<String, String> sequences = new TreeMap<String, String>();

    private int maxLength = 0;
    private int maxWidth = 1024;
    private int ideogramHeight = 20;
    private int ideogramMargin = 25;
    private int xMargin = 40;
    private int yMargin = 40;

    private int getKmerSize(Set<String> kmers) {
        return kmers.iterator().next().length();
    }

    public void initialize() {
        ReferenceSequence rseq;
        while ((rseq = REFERENCE.nextSequence()) != null) {
            String[] names = rseq.getName().split("\\s+");
            String name = names[0];
            String seq = new String(rseq.getBases());

            log.info("loading contig {}...", name);

            if (seq.length() > maxLength) {
                maxLength = seq.length();
            }

            sequences.put(name, seq);
        }
    }

    public void setup() {
        size(3*xMargin + maxWidth, yMargin + sequences.size()*ideogramHeight + sequences.size()*ideogramMargin, PGraphicsPDF.PDF, out.getAbsolutePath());

        noLoop();

        background(Color.WHITE.getRGB());

        //int kmerSize = getKmerSize(KMERS);
        DecimalFormat formatter = new DecimalFormat("#,###");

        int index = 0;
        for (String name : sequences.keySet()) {
            int xpos0 = xMargin;
            int width = (int) (((float) sequences.get(name).length() * maxWidth) / (float) maxLength);
            int ypos0 = yMargin + index*ideogramMargin + index*ideogramHeight;
            int height = ideogramHeight;

            fill(Color.decode("#F0F0F0").getRGB());
            rect(xpos0, ypos0, width, height, 15);

            PFont f = createFont("Helvetica", 15, true);
            textFont(f, 15);
            fill(Color.BLACK.getRGB());
            text(name, xpos0, ypos0 - 4);

            textFont(f, 6);
            fill(Color.BLACK.getRGB());
            text(formatter.format(sequences.get(name).length()) + " bp", xMargin + width + 3, ypos0 + ((float) height / 2.0f) + 2.0f);

            Set<GFF3Record> records = new HashSet<GFF3Record>();

            for (int i = 0; i <= sequences.get(name).length(); i+=WINDOW) {
                SAMRecordIterator it = BAM.queryOverlapping(name, i, i + WINDOW);

                int coverage = 0;
                SAMRecord read;
                while ((read = it.next()) != null) {
                    coverage += read.getMappingQuality() > 0 ? 1 : 0;
                }

                float xpos1 = xMargin + ((float) i*maxWidth/(float) maxLength);
                stroke(Color.decode("#A90641").getRGB());
                strokeWeight(2.0f);



                /*
                String fw = sequences.get(name).substring(i, i + kmerSize);
                String rc = SequenceUtils.reverseComplement(fw);

                if (KMERS.contains(fw) || KMERS.contains(rc)) {
                    float xpos1 = xMargin + ((float) i*maxWidth/(float) maxLength);

                    stroke(Color.decode("#A90641").getRGB());
                    strokeWeight(2.0f);
                    line(xpos1, ypos0, xpos1, ypos0 + height);

                    log.info("kmer: {} contig: {} pos: {}", fw, name, i);

                    Collection<GFF3Record> overlappingRecords = GFF3.getType("gene", GFF.getOverlapping(new Interval(name, i, i)));
                    records.addAll(overlappingRecords);
                }
                */
            }

            for (GFF3Record record : records) {
                float xpos2 = xMargin + ((float) record.getStart()*maxWidth/(float) maxLength);
                float xpos3 = xMargin + ((float) record.getEnd()*maxWidth/(float) maxLength);

                stroke(Color.decode("#0048A9").getRGB());
                strokeWeight(2.5f);
                line(xpos2, ypos0 + height + 4, xpos3, ypos0 + height + 4);

                textFont(f, 3);
                fill(Color.BLACK.getRGB());
                text(record.getAttribute("ID"), xpos3 + 2, ypos0 + height + 5);

                try {
                    text(URLDecoder.decode(record.getAttribute("description"), "UTF-8"), xpos3 + 2, ypos0 + height + 8);
                } catch (UnsupportedEncodingException e) {
                    throw new RuntimeException(e);
                }
            }

            stroke(Color.BLACK.getRGB());
            strokeWeight(2.0f);
            noFill();
            rect(xpos0, ypos0, width, height, 15);

            index++;
        }

        exit();
    }
}
