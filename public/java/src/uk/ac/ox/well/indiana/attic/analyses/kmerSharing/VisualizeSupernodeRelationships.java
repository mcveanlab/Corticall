package uk.ac.ox.well.indiana.attic.analyses.kmerSharing;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import processing.pdf.PGraphicsPDF;
import uk.ac.ox.well.indiana.commands.Sketch;
import uk.ac.ox.well.indiana.utils.alignment.sw.SmithWaterman;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.processing.visualelements.Canvas;
import uk.ac.ox.well.indiana.utils.processing.visualelements.Frame;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.awt.*;
import java.io.File;
import java.util.*;

public class VisualizeSupernodeRelationships extends Sketch {
    @Argument(fullName="reference", shortName="R", doc="Reference sequence")
    public IndexedFastaSequenceFile REFERENCE;

    @Argument(fullName="gff", shortName="gff", doc="GFF file")
    public GFF3 GFF;

    @Argument(fullName="relatedSequences", shortName="rel", doc="Related sequences file")
    public File RELATED_SEQUENCES;

    @Argument(fullName="geneName", shortName="gn", doc="Gene name")
    public String GENE_NAME = "PF3D7_1200600";

    @Output
    public File out;

    private Canvas canvas = new Canvas();

    private class Gene extends Frame {
        private String geneId;
        private GFF3Record gene;
        private Collection<GFF3Record> exons;
        private String seq;
        private HashMap<String, Integer> kmers = new HashMap<String, Integer>();
        private HashMap<String, HashSet<String>> supernodes = new HashMap<String, HashSet<String>>();

        private int xMargin = 10;
        private int yMargin = 30;
        private int exonHeight = 25;
        private int supernodeHeight = 9;

        public Gene(String geneId) {
            this.geneId = geneId;

            gene = GFF.getRecord(geneId);
            exons = GFF.getType("exon", GFF.getOverlapping(gene));

            seq = new String(REFERENCE.getSubsequenceAt(gene.getSeqid(), gene.getStart(), gene.getEnd()).getBases());
        }

        private class SupernodeView implements Comparable<SupernodeView> {
            private String supernode;
            private int pos;
            private HashSet<String> uniqueKmers;
            private SmithWaterman sw;

            public SupernodeView(String supernode, int pos, HashSet<String> uniqueKmers, SmithWaterman sw) {
                this.supernode = supernode;
                this.pos = pos;
                this.uniqueKmers = uniqueKmers;
                this.sw = sw;
            }

            public String getSupernode() { return supernode; }
            public int getPosition() { return pos; }
            public HashSet<String> getUniqueKmers() { return uniqueKmers; }
            public SmithWaterman getSmithWaterman() { return sw; }

            @Override
            public int compareTo(SupernodeView o) {
                if (o.getPosition() == getPosition()) {
                    if (o.getSupernode().length() == getSupernode().length()) {
                        return 0;
                    } else {
                        return o.getSupernode().length() < getSupernode().length() ? -1 : 1;
                    }
                }

                return o.getPosition() < getPosition() ? -1 : 1;
            }
        }

        @Override
        public void draw() {
            int xpos = xMargin;
            int ypos = yMargin;

            fill(Color.BLACK.getRGB());
            text(geneId, xpos, ypos);

            int geneLength = gene.getEnd() - gene.getStart();
            line(xpos, ypos, geneLength + xpos, ypos);

            stroke(Color.BLACK.getRGB());
            fill(Color.WHITE.getRGB());
            for (GFF3Record exon : exons) {
                int exonLength = exon.getEnd() - exon.getStart();
                int exonPos = xpos + exon.getStart() - gene.getStart();

                rect(exonPos, ypos - exonHeight/2, exonLength, exonHeight);
            }

            strokeCap(PROJECT);
            noStroke();
            fill(Color.RED.getRGB(), 50);
            for (String kmer : kmers.keySet()) {
                int pos = xpos + kmers.get(kmer);

                rect(pos, ypos - (exonHeight/2), kmer.length(), exonHeight);
            }

            TreeSet<SupernodeView> supernodeViews = new TreeSet<SupernodeView>();

            for (String fwsupernode : supernodes.keySet()) {
                String rcsupernode = SequenceUtils.reverseComplement(fwsupernode);

                SmithWaterman fsw = new SmithWaterman(seq, fwsupernode);
                SmithWaterman rsw = new SmithWaterman(seq, rcsupernode);

                int pos = 0;
                String goodSupernode = fwsupernode;
                SmithWaterman goodsw = fsw;

                if (fsw.getAlignmentScore() > rsw.getAlignmentScore()) {
                    for (String alignment : fsw.getAlignment()) {
                        if (seq.contains(alignment.replaceAll("-", ""))) {
                            pos = seq.indexOf(alignment.replaceAll("-", ""));
                            goodSupernode = fwsupernode;
                            goodsw = fsw;
                            break;
                        }
                    }
                } else {
                    for (String alignment : rsw.getAlignment()) {
                        if (seq.contains(alignment.replaceAll("-", ""))) {
                            pos = seq.indexOf(alignment.replaceAll("-", ""));
                            goodSupernode = rcsupernode;
                            goodsw = rsw;
                            break;
                        }
                    }
                }

                supernodeViews.add(new SupernodeView(goodSupernode, pos, supernodes.get(fwsupernode), goodsw));
            }

            int supernodeIndex = 0;
            for (SupernodeView sv : supernodeViews) {
                int pos = sv.getPosition();
                int length = sv.getSupernode().length();

                HashSet<Integer> variantPositions = new HashSet<Integer>();

                String alignment1 = sv.getSmithWaterman().getAlignment()[0];
                String alignment2 = sv.getSmithWaterman().getAlignment()[1];

                log.info("{}", alignment1);
                log.info("{}", alignment2);

                int variantPosition = 0;
                for (int i = 0; i < alignment1.length(); i++) {
                    log.info("{} {} {} {}", i, variantPosition, alignment1.charAt(i), alignment2.charAt(i));

                    if (alignment1.charAt(i) == '-') {
                        variantPositions.add(variantPosition);
                    } else if (alignment2.charAt(i) == '-') {
                        variantPositions.add(variantPosition);
                        variantPosition++;
                    } else if (alignment1.charAt(i) != alignment2.charAt(i)) {
                        variantPositions.add(variantPosition);
                        variantPosition++;
                    } else {
                        variantPosition++;
                    }
                }

                log.info("");

                int sxpos = xMargin + pos;
                int sypos = yMargin + exonHeight + (supernodeIndex*supernodeHeight);

                //log.info("{} {}", sv.getSmithWaterman().getAlignment()[0].length(), sv.getSmithWaterman().getAlignment()[1].length());

                noStroke();
                fill(Color.WHITE.getRGB());
                rect(sxpos, sypos, length, (supernodeHeight - 2));

                for (String fwkmer : sv.getUniqueKmers()) {
                    String rckmer = SequenceUtils.reverseComplement(fwkmer);

                    int kmerpos = -1;

                    if (sv.getSupernode().contains(fwkmer)) {
                        kmerpos = sv.getSupernode().indexOf(fwkmer);
                    } else if (sv.getSupernode().contains(rckmer)) {
                        kmerpos = sv.getSupernode().indexOf(rckmer);
                    }

                    if (kmerpos >= 0) {
                        fill(Color.RED.getRGB(), 50);
                        rect(sxpos + kmerpos, sypos, fwkmer.length(), (supernodeHeight - 2));
                    }
                }

                for (int varpos : variantPositions) {
                    stroke(Color.BLUE.getRGB());
                    line(sxpos + varpos, sypos, sxpos + varpos, sypos + supernodeHeight - 2);

                    log.info("variant: {}", varpos);
                }

                supernodeIndex++;
            }
        }

        public void addKmer(String fw) {
            String rc = SequenceUtils.reverseComplement(fw);

            if (seq.contains(fw) || seq.contains(rc)) {
                int pos;

                if (seq.contains(fw)) {
                    pos = seq.indexOf(fw);
                } else {
                    pos = seq.indexOf(rc);
                }

                kmers.put(fw, pos);
            }
        }

        public void addSupernode(String supernode, String fw) {
            String fwsupernode = SequenceUtils.alphanumericallyLowestOrientation(supernode);

            if (!supernodes.containsKey(fwsupernode)) {
                supernodes.put(fwsupernode, new HashSet<String>());
            }

            supernodes.get(fwsupernode).add(fw);
        }

        public int getWidth() {
            return xMargin + (gene.getEnd() - gene.getStart()) + xMargin;
        }

        public int getHeight() {
            return yMargin + (exonHeight + (supernodeHeight * supernodes.size())) + yMargin;
        }
    }

    @Override
    public void initialize() {
        Gene gene = new Gene(GENE_NAME);

        ArrayList<Map<String, String>> relatedSequences = new ArrayList<Map<String, String>>();

        log.info("gene size: {}x{}", gene.getWidth(), gene.getHeight());

        TableReader table = new TableReader(RELATED_SEQUENCES);
        for (Map<String, String> entry : table) {
            if (entry.get("genes").contains(GENE_NAME)) {
                relatedSequences.add(entry);

                gene.addKmer(entry.get("kmer"));
                gene.addSupernode(entry.get("superNode"), entry.get("kmer"));
            }
        }

        log.info("supernodes: {}", gene.supernodes.size());

        canvas.addElement(gene);
    }

    public void setup() {
        size(canvas.getWidth(), canvas.getHeight(), PGraphicsPDF.PDF, out.getAbsolutePath());

        log.info("size: {}x{}", canvas.getWidth(), canvas.getHeight());
    }

    public void draw() {
        canvas.draw();

        exit();
    }
}
