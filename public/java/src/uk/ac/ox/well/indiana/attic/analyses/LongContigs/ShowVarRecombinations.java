package uk.ac.ox.well.indiana.attic.analyses.LongContigs;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import processing.pdf.PGraphicsPDF;
import uk.ac.ox.well.indiana.commands.Sketch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataFrame;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;
import uk.ac.ox.well.indiana.utils.statistics.clustering.HierarchicalClustering;

import java.awt.*;
import java.io.File;
import java.io.PrintStream;
import java.util.*;
import java.util.List;

public class ShowVarRecombinations extends Sketch {
    @Argument(fullName="classifications", shortName="c", doc="Classifications file")
    public LinkedHashMap<String, File> CLASSIFICATIONS;

    @Argument(fullName="reference", shortName="r", doc="Reference FASTA files")
    public TreeMap<String, IndexedFastaSequenceFile> REFERENCES;

    @Argument(fullName="parent", shortName="p", doc="Parent CTX files")
    public TreeMap<String, CortexGraph> PARENTS;

    @Argument(fullName="gff", shortName="g", doc="GFF files")
    public TreeMap<String, GFF3> GFFS;

    @Argument(fullName="colors", shortName="rgb", doc="Color map")
    public HashMap<String, Color> COLORS;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Output
    public File out;

    @Output(fullName="statsOut", shortName="so", doc="Output file for statistics")
    public PrintStream sout;

    private Map<String, Map<String, String>> nameIdContigMap = new TreeMap<String, Map<String, String>>();
    private Map<String, Map<String, Integer>> nameIdPosMap = new HashMap<String, Map<String, Integer>>();
    private Map<String, Map<String, String>> nameIdOriginalNameMap = new HashMap<String, Map<String, String>>();
    private Map<CortexKmer, Color> kmerColorMap = new HashMap<CortexKmer, Color>();
    private Map<CortexKmer, String> kmerNameMap = new HashMap<CortexKmer, String>();

    private int longestContig = 0;
    private int maxContigs = 0;

    private final int marginXRight = 150;
    private final int marginLabel = 160;
    private final int marginX = 10;
    private final int marginY = 10;
    private final int marginTitle = 30;
    private final int marginContig = 5;
    private final int contigHeight = 10;

    public Map<String, Map<String, String>> clusterSequences(List<ReferenceSequence> seqs) {
        DataFrame<String, String, Float> distances = new DataFrame<String, String, Float>(0.0f);
        Map<String, String> seqMap = new HashMap<String, String>();

        for (int i = 0; i < seqs.size(); i++) {
            seqMap.put(seqs.get(i).getName(), new String(seqs.get(i).getBases()));

            for (int j = 0; j < seqs.size(); j++) {
                String seqi = new String(seqs.get(i).getBases());
                Set<CortexKmer> cki = new HashSet<CortexKmer>();
                for (int k = 0; k < seqi.length() - KMER_SIZE; k++) {
                    cki.add(new CortexKmer(seqi.substring(k, k + KMER_SIZE)));
                }

                String seqj = new String(seqs.get(j).getBases());
                Set<CortexKmer> ckj = new HashSet<CortexKmer>();
                for (int k = 0; k < seqj.length() - KMER_SIZE; k++) {
                    ckj.add(new CortexKmer(seqj.substring(k, k + KMER_SIZE)));
                }

                Set<CortexKmer> total = new HashSet<CortexKmer>();
                total.addAll(cki);
                total.addAll(ckj);

                float intersection = 0.0f;
                for (CortexKmer kmer : total) {
                    intersection += (cki.contains(kmer) && ckj.contains(kmer)) ? 1.0f : 0.0f;
                }

                float denom = (float) (cki.size() < ckj.size() ? cki.size() : ckj.size()) + 1.0f;
                float distance = 1.0f - (intersection / denom);

                distances.set(seqs.get(i).getName(), seqs.get(j).getName(), distance);
                distances.set(seqs.get(j).getName(), seqs.get(i).getName(), distance);
            }
        }

        HierarchicalClustering hc = new HierarchicalClustering();
        hc.setMatrix(distances, false);
        //hc.setMembershipLevel(1);
        hc.cluster();
        List<Set<String>> clusters = hc.getClusters();

        Map<String, Map<String, String>> newSeqMap = new TreeMap<String, Map<String, String>>();
        for (int i = 0; i < clusters.size(); i++) {
            String clusterName = String.format("%02d", i);
            newSeqMap.put(clusterName, new HashMap<String, String>());
            nameIdOriginalNameMap.put(clusterName, new HashMap<String, String>());

            for (String m : clusters.get(i)) {
                String[] nameId = m.split("_");
                String contigName = nameId[0];
                String id = nameId[1];

                String contig = seqMap.get(m);

                newSeqMap.get(clusterName).put(id, contig);
                nameIdOriginalNameMap.get(clusterName).put(id, contigName);
            }
        }

        return newSeqMap;
    }

    public void initialize() {
        log.info("Processing classifications...");
        maxContigs = CLASSIFICATIONS.size();

        List<ReferenceSequence> seqs = new ArrayList<ReferenceSequence>();
        List<String> seqStrings = new ArrayList<String>();

        for (String id : CLASSIFICATIONS.keySet()) {
            TableReader tr = new TableReader(CLASSIFICATIONS.get(id));

            int found = 0;
            for (Map<String, String> te : tr) {
                int kmersHB3 = Integer.valueOf(te.get("HB3"));
                int kmers3D7 = Integer.valueOf(te.get("3D7"));

                if (kmersHB3 > 0 && kmers3D7 > 0) {
                    String name = te.get("name");
                    String contig = te.get("contig");

                    seqStrings.add(contig);

                    ReferenceSequence seq = new ReferenceSequence(name + "_" + id, 0, contig.getBytes());
                    seqs.add(seq);

                    found++;
                    if (contig.length() > longestContig) {
                        longestContig = contig.length();
                    }
                }
            }

            log.info("  Found {} chimeric contigs in {}", found, id);
        }

        log.info("  Found {} chimeric contigs total", seqs.size());

        log.info("Finding non-chimeric contigs that overlap chimeras...");
        int count = 0;
        for (String id : CLASSIFICATIONS.keySet()) {
            TableReader tr = new TableReader(CLASSIFICATIONS.get(id));

            for (Map<String, String> te : tr) {
                int kmersHB3 = Integer.valueOf(te.get("HB3"));
                int kmers3D7 = Integer.valueOf(te.get("3D7"));

                if ((kmersHB3 > 0 && kmers3D7 == 0) || (kmersHB3 == 0 && kmers3D7 > 0)) {
                    String name = te.get("name");
                    String contig = te.get("contig");

                    for (String ocontig : seqStrings) {
                        String fw = contig;
                        String rc = SequenceUtils.reverseComplement(contig);

                        if (ocontig.contains(fw) || ocontig.contains(rc)) {
                            ReferenceSequence seq = new ReferenceSequence(name + "_" + id, 0, contig.getBytes());
                            seqs.add(seq);

                            count++;

                            break;
                        }
                    }
                }
            }
        }

        log.info("  Found {} related non-chimeric contigs total", count);

        log.info("Clustering contigs...");
        nameIdContigMap = clusterSequences(seqs);

        log.info("  Found {} clusters", nameIdContigMap.size());

        log.info("Aligning contigs...");
        for (String name : nameIdContigMap.keySet()) {
            String longestContig = "";
            for (String id : nameIdContigMap.get(name).keySet()) {
                String contig = nameIdContigMap.get(name).get(id);

                if (contig.length() > longestContig.length()) {
                    longestContig = contig;
                }
            }

            nameIdPosMap.put(name, new HashMap<String, Integer>());
            for (String id : nameIdContigMap.get(name).keySet()) {
                String fw = nameIdContigMap.get(name).get(id);
                String rc = SequenceUtils.reverseComplement(fw);

                if (longestContig.contains(fw) || longestContig.contains(rc)) {
                    String contig = fw;

                    if (longestContig.contains(rc)) {
                        contig = rc;
                        nameIdContigMap.get(name).put(id, rc);
                    }

                    int offset = longestContig.indexOf(contig);
                    nameIdPosMap.get(name).put(id, offset);
                } else {
                    log.info("  name={} id=={} was not contained in longest contig", name, id);
                    log.info("    l:{}", longestContig);
                    log.info("    f:{}", fw);
                    log.info("    w:{}", rc);

                    nameIdPosMap.get(name).put(id, 0);
                }
            }
        }

        log.info("Painting kmers...");

        for (String refid : REFERENCES.keySet()) {
            log.info("  {} (r={} g={} b={})", refid, COLORS.get(refid).getRed(), COLORS.get(refid).getGreen(), COLORS.get(refid).getBlue());

            IndexedFastaSequenceFile ref = REFERENCES.get(refid);
            CortexGraph cg = PARENTS.get(refid);
            GFF3 gff = GFFS.get(refid);

            for (CortexRecord cr : cg) {
                CortexKmer kmer = cr.getKmer();

                if (!kmerColorMap.containsKey(kmer)) {
                    kmerColorMap.put(kmer, COLORS.get(refid));
                } else {
                    kmerColorMap.put(kmer, Color.LIGHT_GRAY);
                }
            }

            /*
            ReferenceSequence rseq;
            while ((rseq = ref.nextSequence()) != null) {
                String seq = new String(rseq.getBases());

                for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                    CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                    if (!kmerColorMap.containsKey(kmer)) {
                        kmerColorMap.put(kmer, COLORS.get(refid));
                    } else {
                        kmerColorMap.put(kmer, Color.LIGHT_GRAY);
                    }
                }
            }
            */

            for (GFF3Record gr : gff) {
                if (gr.getType().equals("gene")) {
                    String geneseq = new String(ref.getSubsequenceAt(gr.getSeqid(), gr.getStart(), gr.getEnd()).getBases());

                    for (int i = 0; i <= geneseq.length() - KMER_SIZE; i++) {
                        CortexKmer kmer = new CortexKmer(geneseq.substring(i, i + KMER_SIZE));

                        if (!kmerNameMap.containsKey(kmer)) {
                            kmerNameMap.put(kmer, gr.getAttribute("ID") + " (" + gr.getSeqid() + ":" + gr.getStart() + "-" + gr.getEnd() + ")");
                        } else {
                            kmerNameMap.put(kmer, "");
                        }
                    }
                }
            }
        }
    }

    public void setup() {
        size(marginLabel + 2*marginX + longestContig + marginXRight, marginTitle + 2*marginY + (marginContig + contigHeight)*(maxContigs + 10), PGraphicsPDF.PDF, out.getAbsolutePath());

        background(Color.WHITE.getRGB());

        PGraphicsPDF pdf = (PGraphicsPDF) g;

        Map<String, String> labelMapping = new HashMap<String, String>();
        labelMapping.put("sn", "supernode");
        labelMapping.put("se", "single-end");
        labelMapping.put("pe1", "paired-end (one-way)");
        labelMapping.put("pe2", "paired-end (two-way)");

        List<String> order = new ArrayList<String>();
        order.add("pe2");
        order.add("pe1");
        order.add("se");
        order.add("sn");

        Map<String, Integer> longestContigMap = new HashMap<String, Integer>();
        longestContigMap.put("sn", 0);
        longestContigMap.put("se", 0);
        longestContigMap.put("pe1", 0);
        longestContigMap.put("pe2", 0);
        int equalPe1Pe2 = 0;

        log.info("Making plots...");

        TableWriter tw = new TableWriter(sout);

        int index = 0;
        for (String name : nameIdContigMap.keySet()) {
            if (nameIdContigMap.get(name).size() > 1) {
                log.info("  page {}", index);
                if (index > 0) { pdf.nextPage(); }
                index++;

                background(Color.WHITE.getRGB());

                fill(Color.BLACK.getRGB());
                textAlign(LEFT, CENTER);
                textSize(15);
                text("Cluster " + name, marginX + 30, marginTitle / 2);

                int cxpos = marginLabel + marginX;
                int cypos = marginTitle + marginY;
                String longestContig = "";
                for (String id : nameIdContigMap.get(name).keySet()) {
                    String contig = nameIdContigMap.get(name).get(id);

                    if (contig.length() > longestContig.length()) {
                        longestContig = contig;
                    }
                }

                textSize(11);
                textAlign(LEFT, TOP);
                fill(Color.BLACK.getRGB());
                text("genes", marginX + 30, cypos + 5*(marginContig + contigHeight));
                //int ypos = marginTitle + marginY + i*(marginContig + contigHeight);

                String prevName = null;
                int prevYFactor = 0;
                int prevYPlace = 0;
                for (int p = 0; p < longestContig.length(); p++) {
                    char base = longestContig.charAt(p);

                    Color baseColor;

                    switch (base) {
                        case 'A': baseColor = Color.GREEN;  break;
                        case 'C': baseColor = Color.BLUE;   break;
                        case 'G': baseColor = Color.ORANGE; break;
                        case 'T': baseColor = Color.RED;    break;
                        default:  baseColor = Color.BLACK;  break;
                    }

                    stroke(baseColor.getRGB(), 255.0f);
                    strokeCap(PROJECT);
                    line(cxpos + p, cypos + 1, cxpos + p, cypos + contigHeight - 1);

                    if (p < longestContig.length() - KMER_SIZE) {
                        CortexKmer kmer = new CortexKmer(longestContig.substring(p, p + KMER_SIZE));
                        if (kmerNameMap.containsKey(kmer)) {
                            String currentName = kmerNameMap.get(kmer);

                            if (!currentName.isEmpty() && (prevName == null || !prevName.equals(currentName))) {
                                if (p - prevYPlace > 100) {
                                    prevYFactor = 0;
                                }

                                textAlign(LEFT, BOTTOM);
                                textSize(6);
                                text(currentName, cxpos + p, cypos + (6 + prevYFactor)*(marginContig + contigHeight));

                                log.info("    {} {} {}", currentName, prevYFactor, prevYPlace);

                                prevYFactor++;
                                prevYPlace = p;

                                prevName = currentName;
                            }
                        }
                    }
                }

                String longestId = "sn";
                int longestLength = 0;
                int pe1Length = 0;
                int pe2Length = 0;

                Map<String, String> statsEntry = new LinkedHashMap<String, String>();
                statsEntry.put("sn", "0");
                statsEntry.put("se", "0");
                statsEntry.put("pe1", "0");
                statsEntry.put("pe2", "0");

                int i = 1;
                for (String id : order) {
                    if (nameIdContigMap.get(name).containsKey(id)) {
                        String contig = nameIdContigMap.get(name).get(id);

                        if (contig.length() >= longestLength) {
                            longestId = id;
                            longestLength = contig.length();
                        }

                        if (id.equals("pe1")) { pe1Length = contig.length(); }
                        if (id.equals("pe2")) { pe2Length = contig.length(); }

                        statsEntry.put(id, String.valueOf(contig.length()));

                        int offset = 0;
                        if (nameIdPosMap.containsKey(name) && nameIdPosMap.get(name).containsKey(id)) {
                            offset = nameIdPosMap.get(name).get(id);
                        }

                        int xpos = marginLabel + marginX + offset;
                        int ypos = marginTitle + marginY + i*(marginContig + contigHeight);

                        textSize(11);
                        textAlign(LEFT, TOP);
                        fill(Color.BLACK.getRGB());
                        text(labelMapping.get(id), marginX + 30, ypos - 2);

                        textSize(10);
                        textAlign(LEFT, TOP);
                        fill(Color.BLACK.getRGB());
                        text(contig.length() + " bp (" + nameIdOriginalNameMap.get(name).get(id) + ")", xpos + contig.length() + 3, ypos - 2);

                        for (int k = 0; k <= contig.length() - KMER_SIZE; k++) {
                            CortexKmer kmer = new CortexKmer(contig.substring(k, k + KMER_SIZE));

                            if (kmerColorMap.containsKey(kmer) && kmerColorMap.get(kmer).equals(Color.LIGHT_GRAY)) {
                                Color color = kmerColorMap.get(kmer);

                                stroke(color.getRGB(), 255.0f);
                                strokeCap(PROJECT);

                                for (int q = 0; q < KMER_SIZE; q++) {
                                    line(xpos + k + q, ypos + 1, xpos + k + q, ypos + contigHeight - 1);
                                }
                            }
                        }

                        for (int k = 0; k <= contig.length() - KMER_SIZE; k++) {
                            CortexKmer kmer = new CortexKmer(contig.substring(k, k + KMER_SIZE));

                            if (kmerColorMap.containsKey(kmer) && !kmerColorMap.get(kmer).equals(Color.LIGHT_GRAY)) {
                                Color color = kmerColorMap.get(kmer);

                                stroke(color.getRGB(), 255.0f);
                                strokeCap(PROJECT);

                                for (int q = 0; q < KMER_SIZE; q++) {
                                    line(xpos + k + q, ypos + 1, xpos + k + q, ypos + contigHeight - 1);
                                }
                            }
                        }

                        stroke(Color.BLACK.getRGB());
                        noFill();
                        rect(xpos - 1, ypos, contig.length() + 1, contigHeight);

                        i++;
                    }
                }

                tw.addEntry(statsEntry);

                longestContigMap.put(longestId, longestContigMap.get(longestId) + 1);
                if (pe1Length == pe2Length) {
                    equalPe1Pe2++;
                }
            }
        }

        log.info("Stats:");
        log.info("     sn: {}", longestContigMap.get("sn"));
        log.info("     se: {}", longestContigMap.get("se"));
        log.info("    pe1: {}", longestContigMap.get("pe1"));
        log.info("    pe2: {}", longestContigMap.get("pe2"));
        log.info("pe1=pe2: {}", equalPe1Pe2);
    }

    public void draw() {
        exit();
    }
}
