package uk.ac.ox.well.indiana.attic.analyses.reconstruction;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import processing.pdf.PGraphicsPDF;
import uk.ac.ox.well.indiana.commands.Sketch;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.awt.*;
import java.io.File;
import java.util.*;
import java.util.List;

public class ViewVarRecombs extends Sketch {
    @Argument(fullName="graph", shortName="g", doc="Trio graph")
    public CortexGraph GRAPH;

    @Argument(fullName="ref", shortName="r", doc="Reference genome")
    public File REF;

    @Argument(fullName="proteinDomains", shortName="pd", doc="Protein domains")
    public File PROTEIN_DOMAINS;

    @Argument(fullName="genes", shortName="gn", doc="Genes to process")
    public FastaSequenceFile GENES;

    @Argument(fullName="geneOrder", shortName="go", doc="Gene order", required=false)
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

    private Map<String, Color> geneColor = new HashMap<String, Color>();
    private Map<String, String> geneRecords = new HashMap<String, String>();
    private Map<String, GeneView> genes = new TreeMap<String, GeneView>();
    private Map<CortexKmer, Color> kmers = new HashMap<CortexKmer, Color>();

    private KmerLookup kl;

    private boolean isUnique(CortexKmer kmer) {
        String fw = kmer.getKmerAsString();
        String rc = SequenceUtils.reverseComplement(fw);

        return (kl.findKmer(fw).size() == 1 && kl.findKmer(rc).size() == 0) || (kl.findKmer(fw).size() == 0 && kl.findKmer(rc).size() == 1);
    }

    private class GeneView {
        public String geneName;
        public String geneCds;
//        public GFF3Record geneRecord;
//        public Collection<GFF3Record> records;
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

            int proteinWidth = geneCds.length() / 3;

            // Print a line denoting the length of the protein
            strokeWeight(1);
            strokeCap(SQUARE);
            stroke(Color.BLACK.getRGB());
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

            for (int j = 0, codonPos = 0, aaPos = 0; j < geneCds.length() - GRAPH.getKmerSize(); j++, codonPos++) {
                if (codonPos >= 3) {
                    aaPos++;
                    codonPos = 0;
                }

                if (!aaToGenomic.containsKey(aaPos)) {
                    aaToGenomic.put(aaPos, j);
                }
            }

            int kmersRecovered = 0, kmersUnrecovered = 0;

            // Draw kmers
            for (int aaPos : aaToGenomic.keySet()) {
                int gPos = aaToGenomic.get(aaPos);

                CortexKmer kmer = new CortexKmer(geneCds.substring(gPos, gPos + GRAPH.getKmerSize()));
                CortexRecord cr = GRAPH.findRecord(kmer);

                if (isUnique(kmer)) {
                    if (cr != null && cr.getCoverage(0) > 0 && kmers.containsKey(kmer)) {
                        stroke(kmers.get(kmer).getRGB(), 150.0f);
                        strokeWeight(1);
                        strokeCap(SQUARE);
                        line(x + labelMargin + aaPos, y, x + labelMargin + aaPos, y + height);

                        kmersRecovered++;
                    } else {
                        stroke(Color.LIGHT_GRAY.getRGB(), 60.0f);
                        strokeWeight(1);
                        strokeCap(SQUARE);
                        line(x + labelMargin + aaPos, y, x + labelMargin + aaPos, y + height);

                        kmersUnrecovered++;
                    }
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

                    strokeWeight(1);
                    stroke(Color.BLACK.getRGB());
                    line(x + labelMargin + domainOffset + domainWidth - 8, y + (height/2), x + labelMargin + domainOffset + domainWidth, y + height/2);
                } else if ("ATS".equalsIgnoreCase(pd.domainClass)) {
                    stroke(Color.WHITE.getRGB());
                    for (int i = 0; i < 20; i++) {
                        quad(x + labelMargin + domainOffset, y, x + labelMargin + domainOffset + domainWidth, y, x + labelMargin + domainOffset + domainWidth - i, y + height, x + labelMargin + domainOffset + i, y + height);
                    }

                    strokeWeight(1);
                    stroke(Color.BLACK.getRGB());
                    line(x + labelMargin + domainOffset - 5, y + (height/2), x + labelMargin + domainOffset + 10, y + height/2);
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

            // Print gene info
            float pctRecovered = 100.0f * ((float) kmersRecovered) / ((float) (kmersRecovered + kmersUnrecovered));
            //String caption = geneRecord.getSeqid() + ", " + geneRecord.getAttribute("position") + ", " + geneRecord.getAttribute("class") + ", " + String.format("%.2f%%", pctRecovered);
            String caption = geneName;
            textSize(18);
            textAlign(LEFT, CENTER);
            text(caption, x, y + 70);
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
            rect(x + labelMargin + domainOffset, y, domainWidth, height);

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
            String newid = te.get("oldid");

            if (!proteinDomains.containsKey(newid)) {
                proteinDomains.put(newid, new ArrayList<ProteinDomain>());
            }

            ProteinDomain pd = new ProteinDomain();
            //pd.domainClass = te.get("class");

            if (te.get("domain").startsWith("NTS")) {
                pd.domainClass = "NTS";
            } else if (te.get("domain").startsWith("ATS")) {
                pd.domainClass = "ATS";
            } else if (te.get("domain").startsWith("CIDR")) {
                pd.domainClass = "CIDR";
            } else if (te.get("domain").startsWith("DBL")) {
                pd.domainClass = "DBL";
            } else {
                pd.domainClass = "UKN";
            }

            pd.domainName = te.get("domain");
            pd.aaStart = Integer.valueOf(te.get("aa_start"));
            pd.aaEnd = Integer.valueOf(te.get("aa_end"));

            proteinDomains.get(newid).add(pd);
        }

        return proteinDomains;
    }

    private ArrayList<String> getGeneNames() {
        if (GENE_ORDER == null) {
            Set<String> gos = new TreeSet<String>();

            TableReader tr = new TableReader(PROTEIN_DOMAINS);
            for (Map<String, String> te : tr) {
                gos.add(te.get("oldid"));
            }

            return new ArrayList<String>(gos);
        }

        return GENE_ORDER;
    }

    public void initialize() {
        kl = new KmerLookup(REF, GRAPH.getKmerSize());

        // Load gene order
        ArrayList<String> geneNames = getGeneNames();

        Color[] colors = generateColors(geneNames.size());
        for (int i = 0; i < geneNames.size(); i++) {
            geneColor.put(geneNames.get(i), colors[i]);
        }

        // Load gene records
        Map<String, Collection<ProteinDomain>> proteinDomains = loadProteinDomains();
        ReferenceSequence rseq;
        while ((rseq = GENES.nextSequence()) != null) {
            String id = rseq.getName();

            if (geneColor.containsKey(id)) {
                String seq = new String(rseq.getBases());

                geneRecords.put(id, seq);

                int recovery = 0;
                for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                    CortexKmer kmer = new CortexKmer(seq.substring(i, i + GRAPH.getKmerSize()));
                    CortexRecord cr = GRAPH.findRecord(kmer);

                    if (cr != null && cr.getCoverage(0) > 0) {
                        recovery++;
                    }

                    if (!kmers.containsKey(kmer) && isUnique(kmer)) {
                        kmers.put(kmer, geneColor.get(id));
                    }
                }

                log.info("{}: {}/{}", id, recovery, seq.length() - GRAPH.getKmerSize());
            }
        }

        // Display genes
        int longestWidth = 0;
        for (String geneName : geneNames) {
            if (geneRecords.containsKey(geneName)) {
                String seq = geneRecords.get(geneName);

                GeneView gv = new GeneView();
                gv.geneName = geneName;
                gv.geneCds = seq;
                gv.proteinDomains = proteinDomains.get(geneName);

                gv.width = Math.round(((float) (seq.length())) / 3.0f);

                if (gv.width > longestWidth) { longestWidth = gv.width; }

                genes.put(geneName, gv);
            }
        }

        windowWidth = 2*windowXMargin + geneLabelMargin + longestWidth + 2*legendMargin;
        windowHeight = 2*windowYMargin + geneMargin*geneRecords.size();
    }

    public void setup() {
        size(windowWidth, windowHeight, PGraphicsPDF.PDF, out.getAbsolutePath());

        background(Color.WHITE.getRGB());
    }

    public void draw() {
        int geneIndex = 0;
        for (String gene : getGeneNames()) {
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
