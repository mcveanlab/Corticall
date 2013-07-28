package uk.ac.ox.well.indiana.analyses.reconstruction;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTreeMap;
import processing.pdf.PGraphicsPDF;
import uk.ac.ox.well.indiana.sketches.Sketch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

import java.awt.*;
import java.io.File;
import java.util.*;
import java.util.List;

public class ViewKmerPlot extends Sketch {
    @Argument(fullName="reference", shortName="R", doc="Reference genome")
    public IndexedFastaSequenceFile FASTA;

    @Argument(fullName="gff", shortName="gff", doc="GFF of regions of interest")
    public GFF3 GFF;

    @Argument(fullName="kmerReferencePanel", shortName="krp", doc="Kmer reference panel to color")
    public File KMER_REFERENCE_PANEL;

    @Argument(fullName="contigsTable", shortName="ct", doc="Contigs table")
    public File CONTIGS_TABLE;

    @Argument(fullName="proteinDomains", shortName="pd", doc="Protein domains")
    public File PROTEIN_DOMAINS;

    @Argument(fullName="geneOrder", shortName="go", doc="Gene order")
    public ArrayList<String> GENE_ORDER;

    @Output
    public File out;

    private int windowXMargin = 50;
    private int windowYMargin = 50;
    private int geneMargin = 150;
    private int geneHeight = 50;
    private int geneLabelMargin = 450;
    private int legendMargin = 300;
    private int legendElementHeight = 30;

    private int windowWidth;
    private int windowHeight;
    private int kmerSize;

    private Map<String, Color> geneColor = new HashMap<String, Color>();
    private Map<String, GFF3Record> geneRecords = new HashMap<String, GFF3Record>();
    private Map<String, GeneView> genes = new TreeMap<String, GeneView>();
    private Map<CortexKmer, Color> kmers;
    private Map<CortexKmer, Color> kmerReferencePanel;

    private class GeneView {
        public String geneName;
        public GFF3Record geneRecord;
        public Collection<GFF3Record> records;
        public Collection<ProteinDomain> proteinDomains;
        public IntervalTreeMap<ProteinDomain> proteinDomainIntervals;
        public int width;

        public void draw(int x, int y, int labelMargin, int height) {
            proteinDomainIntervals = new IntervalTreeMap<ProteinDomain>();

            // Print gene name
            textSize(50);
            textAlign(LEFT, CENTER);
            fill(Color.BLACK.getRGB());
            text(geneName, x, y + 20);

            int proteinWidth = 0;
            for (ProteinDomain pd : proteinDomains) {
                if (pd.aaEnd > proteinWidth) {
                    proteinWidth = pd.aaEnd;
                }
            }

            // Print a line denoting the length of the protein
            strokeWeight(1);
            strokeCap(SQUARE);
            fill(Color.GRAY.getRGB());
            line(x + labelMargin, y + (height/2), x + labelMargin + proteinWidth, y + (height/2));

            // Print the length of the protein
            textSize(25);
            textAlign(LEFT, CENTER);
            fill(Color.BLACK.getRGB());
            text(proteinWidth + "aa", x + labelMargin + proteinWidth + 10, y + (height/2) - 3);

            // For each domain, draw an opaque element denoting its boundaries
            for (ProteinDomain pd : proteinDomains) {
                fill(Color.WHITE.getRGB());
                stroke(Color.BLACK.getRGB());
                strokeWeight(5);

                int domainOffset = pd.aaStart;
                int domainWidth = pd.aaEnd - pd.aaStart;

                if ("DBL".equalsIgnoreCase(pd.domainClass)) {
                    rect(x + labelMargin + domainOffset, y, domainWidth, height, 50);
                } else if ("CIDR".equalsIgnoreCase(pd.domainClass) || "CIDRpam".equalsIgnoreCase(pd.domainClass)) {
                    rect(x + labelMargin + domainOffset, y, domainWidth, height);
                } else if ("NTS".equalsIgnoreCase(pd.domainClass) || "NTSpam".equalsIgnoreCase(pd.domainClass)) {
                    quad(x + labelMargin + domainOffset + 20, y, x + labelMargin + domainOffset + domainWidth - 20, y, x + labelMargin + domainOffset + domainWidth, y + height, x + labelMargin + domainOffset, y + height);
                } else if ("ATS".equalsIgnoreCase(pd.domainClass)) {
                    quad(x + labelMargin + domainOffset, y, x + labelMargin + domainOffset + domainWidth, y, x + labelMargin + domainOffset + domainWidth - 20, y + height, x + labelMargin + domainOffset + 20, y + height);
                }
            }

            // Determine mapping from kmer's genomic coordinates to protein coordinates
            Map<Integer, Integer> aaToGenomic = new TreeMap<Integer, Integer>();

            Map<Interval, GFF3Record> sortedRecords = new TreeMap<Interval, GFF3Record>();
            for (GFF3Record record : GFF.getType("exon", records)) {
                sortedRecords.put(record.getInterval(), record);
            }

            List<GFF3Record> recs = new ArrayList<GFF3Record>(sortedRecords.values());
            if (geneRecord.getStrand() == GFF3Record.Strand.POSITIVE) {
                int aaPos = 1;
                for (int i = 0; i < recs.size(); i++) {
                    GFF3Record exon = recs.get(i);
                    int exonLength = exon.getEnd() - exon.getStart();

                    for (int j = 0, codonPos = 0; j < exonLength; j++, codonPos++) {
                        if (codonPos >= 3) {
                            aaPos++;
                            codonPos = 0;
                        }
                        int gPos = exon.getStart() + j;

                        if (!aaToGenomic.containsKey(aaPos)) {
                            aaToGenomic.put(aaPos, gPos);
                        }
                    }
                }
            } else {
                int aaPos = 1;
                for (int i = recs.size() - 1; i >= 0; i--) {
                    GFF3Record exon = recs.get(i);
                    int exonLength = exon.getEnd() - exon.getStart();

                    for (int j = 0, codonPos = 0; j < exonLength; j++, codonPos++) {
                        if (codonPos >= 3) {
                            aaPos++;
                            codonPos = 0;
                        }
                        int gPos = exon.getEnd() - j;

                        if (!aaToGenomic.containsKey(aaPos)) {
                            aaToGenomic.put(aaPos, gPos);
                        }
                    }
                }
            }

            // Draw kmers
            for (int aaPos : aaToGenomic.keySet()) {
                int gPos = aaToGenomic.get(aaPos);

                CortexKmer kmer;

                if (geneRecord.getStrand() == GFF3Record.Strand.POSITIVE) {
                    kmer = new CortexKmer(FASTA.getSubsequenceAt(geneRecord.getSeqid(), gPos, gPos + kmerSize - 1).getBases());
                } else {
                    kmer = new CortexKmer(FASTA.getSubsequenceAt(geneRecord.getSeqid(), gPos - kmerSize + 1, gPos).getBases());
                }

                if (kmers.containsKey(kmer)) {
                    stroke(kmers.get(kmer).getRGB(), 150.0f);
                    strokeWeight(1);
                    strokeCap(SQUARE);
                    line(x + labelMargin + aaPos, y, x + labelMargin + aaPos, y + height);
                } else if (kmerReferencePanel.containsKey(kmer)) {
                    stroke(kmerReferencePanel.get(kmer).getRGB(), 30.0f);
                    strokeWeight(1);
                    strokeCap(SQUARE);
                    line(x + labelMargin + aaPos, y, x + labelMargin + aaPos, y + height);
                }
            }

            // Wipe away parts of kmer lines that fall outside the protein domain bounding boxes
            for (ProteinDomain pd : proteinDomains) {
                noFill();
                stroke(Color.WHITE.getRGB());
                strokeWeight(5);

                int domainOffset = pd.aaStart;
                int domainWidth = pd.aaEnd - pd.aaStart;

                if ("DBL".equalsIgnoreCase(pd.domainClass)) {
                    for (int i = 0; i < 50; i++) {
                        rect(x + labelMargin + domainOffset, y, domainWidth, height, i);
                    }
                } else if ("CIDR".equalsIgnoreCase(pd.domainClass) || "CIDRpam".equalsIgnoreCase(pd.domainClass)) {
                } else if ("NTS".equalsIgnoreCase(pd.domainClass) || "NTSpam".equalsIgnoreCase(pd.domainClass)) {
                    for (int i = 0; i < 20; i++) {
                        quad(x + labelMargin + domainOffset + i, y, x + labelMargin + domainOffset + domainWidth - i, y, x + labelMargin + domainOffset + domainWidth, y + height, x + labelMargin + domainOffset, y + height);
                    }
                } else if ("ATS".equalsIgnoreCase(pd.domainClass)) {
                    stroke(Color.WHITE.getRGB());
                    for (int i = 0; i < 20; i++) {
                        quad(x + labelMargin + domainOffset, y, x + labelMargin + domainOffset + domainWidth, y, x + labelMargin + domainOffset + domainWidth - i, y + height, x + labelMargin + domainOffset + i, y + height);
                    }
                }
            }

            // Redraw protein domain boundaries
            for (ProteinDomain pd : proteinDomains) {
                noFill();
                stroke(Color.BLACK.getRGB());
                strokeWeight(5);

                int domainOffset = pd.aaStart;
                int domainWidth = pd.aaEnd - pd.aaStart;

                if ("DBL".equalsIgnoreCase(pd.domainClass)) {
                    rect(x + labelMargin + domainOffset, y, domainWidth, height, 50);
                } else if ("CIDR".equalsIgnoreCase(pd.domainClass) || "CIDRpam".equalsIgnoreCase(pd.domainClass)) {
                    rect(x + labelMargin + domainOffset, y, domainWidth, height);
                } else if ("NTS".equalsIgnoreCase(pd.domainClass) || "NTSpam".equalsIgnoreCase(pd.domainClass)) {
                    quad(x + labelMargin + domainOffset + 20, y, x + labelMargin + domainOffset + domainWidth - 20, y, x + labelMargin + domainOffset + domainWidth, y + height, x + labelMargin + domainOffset, y + height);
                } else if ("ATS".equalsIgnoreCase(pd.domainClass)) {
                    quad(x + labelMargin + domainOffset, y, x + labelMargin + domainOffset + domainWidth, y, x + labelMargin + domainOffset + domainWidth - 20, y + height, x + labelMargin + domainOffset + 20, y + height);
                }

                textSize(25);
                textAlign(CENTER, CENTER);
                fill(Color.BLACK.getRGB());
                text(pd.domainName, x + labelMargin + domainOffset + (domainWidth/2), y + height + 15);
            }
        }
    }

    private class Legend {
        private int domainWidth = 100;
        private int domainOffset = 100;
        private int labelMargin = 50;

        public void draw(int x, int y, int height) {
            textSize(50);

            fill(Color.BLACK.getRGB());
            text("NTS", x, y + 10);
            noFill();
            quad(x + labelMargin + domainOffset + 20, y, x + labelMargin + domainOffset + domainWidth - 20, y, x + labelMargin + domainOffset + domainWidth, y + height, x + labelMargin + domainOffset, y + height);

            y += 50;
            fill(Color.BLACK.getRGB());
            text("DBL", x, y + 10);
            noFill();
            rect(x + labelMargin + domainOffset, y, domainWidth, height, 50);

            y += 50;
            fill(Color.BLACK.getRGB());
            text("CIDR", x, y + 10);
            noFill();
            rect(x + labelMargin + domainOffset, y - 5, domainWidth, height);

            y += 50;
            fill(Color.BLACK.getRGB());
            text("ATS", x, y + 10);
            noFill();
            quad(x + labelMargin + domainOffset, y, x + labelMargin + domainOffset + domainWidth, y, x + labelMargin + domainOffset + domainWidth - 20, y + height, x + labelMargin + domainOffset + 20, y + height);
        }
    }

    private class ProteinDomain {
        public String domainClass;
        public String domainName;
        public int aaStart;
        public int aaEnd;

        public String toString() {
            return String.format("domainClass=%s domainName=%s aaStart=%d aaEnd=%d", domainClass, domainName, aaStart, aaEnd);
        }
    }

    private Color[] generateColors(int n) {
        Color[] cols = new Color[n];

        for(int i = 0; i < n; i++) {
            cols[i] = Color.getHSBColor((float) i / (float) n, 0.85f, 1.0f);
        }

        return cols;
    }

    private Map<String, Collection<ProteinDomain>> loadProteinDomains() {
        Map<String, Collection<ProteinDomain>> proteinDomains = new HashMap<String, Collection<ProteinDomain>>();

        TableReader tr = new TableReader(PROTEIN_DOMAINS);

        for (Map<String, String> te : tr) {
            String newid = te.get("newid");

            if (!proteinDomains.containsKey(newid)) {
                proteinDomains.put(newid, new ArrayList<ProteinDomain>());
            }

            ProteinDomain pd = new ProteinDomain();
            pd.domainClass = te.get("class");
            pd.domainName = te.get("domain");
            pd.aaStart = Integer.valueOf(te.get("aa_start"));
            pd.aaEnd = Integer.valueOf(te.get("aa_end"));

            proteinDomains.get(newid).add(pd);
        }

        return proteinDomains;
    }

    private Map<CortexKmer, Color> loadContigsTable(File table) {
        Map<CortexKmer, Color> krp = new HashMap<CortexKmer, Color>();

        TableReader tr = new TableReader(table);

        for (Map<String, String> entry : tr) {
            if (entry.containsKey("kmer")) {
                krp.put(new CortexKmer(entry.get("kmer")), Color.LIGHT_GRAY);
            } else {
                String[] kmers = entry.get("kmers").split(",");

                for (String kmer : kmers) {
                    krp.put(new CortexKmer(kmer), Color.LIGHT_GRAY);
                }
            }
        }

        return krp;
    }

    public void initialize() {
        // Load gene order
        Color[] colors = generateColors(GENE_ORDER.size());
        for (int i = 0; i < GENE_ORDER.size(); i++) {
            geneColor.put(GENE_ORDER.get(i), colors[i]);
        }

        // Load kmers
        kmerReferencePanel = loadContigsTable(KMER_REFERENCE_PANEL);
        kmers = loadContigsTable(CONTIGS_TABLE);
        kmerSize = kmers.keySet().iterator().next().length();

        // Load gene records
        Map<String, Collection<ProteinDomain>> proteinDomains = loadProteinDomains();
        for (GFF3Record gr : GFF) {
            String id = gr.getAttribute("ID");

            if (geneColor.containsKey(id)) {
                geneRecords.put(id, gr);

                String seq = new String(FASTA.getSubsequenceAt(gr.getSeqid(), gr.getStart(), gr.getEnd()).getBases());
                for (int i = 0; i <= seq.length() - kmerSize; i++) {
                    CortexKmer kmer = new CortexKmer(seq.substring(i, i + kmerSize));

                    if (kmers.containsKey(kmer)) {
                        kmers.put(kmer, geneColor.get(id));
                    }
                }
            }
        }

        // Display genes
        int longestWidth = 0;
        for (String geneName : GENE_ORDER) {
            if (geneRecords.containsKey(geneName)) {
                GFF3Record gr = geneRecords.get(geneName);

                GeneView gv = new GeneView();
                gv.geneName = geneName;
                gv.geneRecord = gr;
                gv.records = GFF.getOverlapping(gr);
                gv.proteinDomains = proteinDomains.get(gr.getAttribute("ID"));
                gv.width = 0;

                for (ProteinDomain pd : gv.proteinDomains) {
                    if (pd.aaEnd > gv.width) {
                        gv.width = pd.aaEnd;
                    }
                }

                if (gv.width > longestWidth) { longestWidth = gv.width; }

                genes.put(geneName, gv);
            }
        }

        windowWidth = 2*windowXMargin + geneLabelMargin + longestWidth + legendMargin;
        windowHeight = 2*windowYMargin + geneMargin*geneRecords.size();
    }

    public void setup() {
        size(windowWidth, windowHeight, PGraphicsPDF.PDF, out.getAbsolutePath());

        background(Color.WHITE.getRGB());
    }

    public void draw() {
        int geneIndex = 0;
        for (String gene : GENE_ORDER) {
            if (genes.containsKey(gene)) {
                GeneView gv = genes.get(gene);

                int x = windowXMargin;
                int y = windowYMargin + geneIndex*geneMargin;

                gv.draw(x, y, geneLabelMargin, geneHeight);

                geneIndex++;
            }
        }

        Legend legend = new Legend();
        legend.draw(windowWidth - legendMargin, windowYMargin, legendElementHeight);

        exit();
    }
}
