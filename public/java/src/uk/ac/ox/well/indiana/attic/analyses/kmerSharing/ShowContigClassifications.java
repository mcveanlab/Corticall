package uk.ac.ox.well.indiana.attic.analyses.kmerSharing;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.RExecutor;
import processing.pdf.PGraphicsPDF;
import uk.ac.ox.well.indiana.commands.Sketch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataFrame;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

public class ShowContigClassifications extends Sketch {
    @Argument(fullName="classifications", shortName="c", doc="Classifications file")
    public LinkedHashMap<String, File> CLASSIFICATIONS;

    @Argument(fullName="reference", shortName="r", doc="Reference FASTA files")
    public TreeMap<String, IndexedFastaSequenceFile> REFERENCES;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Output
    public File out;

    @Output(fullName="mo", shortName="mo", doc="Distance matrix")
    public File m;

    @Output(fullName="go", shortName="go", doc="Groupings")
    public File g1;

    private Map<String, Map<String, String>> nameIdContigMap = new HashMap<String, Map<String, String>>();
    private Map<String, Map<String, Integer>> nameIdPosMap = new HashMap<String, Map<String, Integer>>();
    private Map<CortexKmer, Color> kmerColorMap = new HashMap<CortexKmer, Color>();

    private final int marginXRight = 50;
    private final int marginLabel = 160;
    private final int marginX = 10;
    private final int marginY = 10;
    private final int marginTitle = 30;
    private final int marginContig = 5;
    private final int contigHeight = 10;
    private int longestContig = 0;
    private int maxContigs = 0;

    private final String clusterGenesScript = "R/clustergenes.R";

    public Map<String, Map<String, String>> clusterSequences(Map<String, Map<String, String>> seqMap) {
        List<ReferenceSequence> seqs = new ArrayList<ReferenceSequence>();

        for (String name : seqMap.keySet()) {
            for (String id : seqMap.get(name).keySet()) {
                String contig = seqMap.get(name).get(id);

                ReferenceSequence rseq = new ReferenceSequence(name + "_" + id, 0, contig.getBytes());

                log.info("{}", rseq.getName());

                seqs.add(rseq);
            }
        }

        DataFrame<String, String, Float> distances = new DataFrame<String, String, Float>(0.0f);

        for (int i = 0; i < seqs.size(); i++) {
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

                float denom = (float) (cki.size() < ckj.size() ? cki.size() : ckj.size());
                float distance = 1.0f - (intersection / denom);

                distances.set(seqs.get(i).getName(), seqs.get(j).getName(), distance);
                distances.set(seqs.get(j).getName(), seqs.get(i).getName(), distance);
            }
        }

        Map<String, Map<String, String>> newSeqMap = new HashMap<String, Map<String, String>>();

        try {
            PrintStream mout = new PrintStream(m);
            mout.print(distances);

            RExecutor.executeFromClasspath(clusterGenesScript, m.getAbsolutePath(), g1.getAbsolutePath(), "-1");

            TableReader tr = new TableReader(g1);
            for (Map<String, String> te : tr) {
                log.info("te:{}", te);

                String[] members = te.get("groupMembers").split(",");

                String newname = "group" + te.get("groupName");

                if (!newSeqMap.containsKey(newname)) {
                    newSeqMap.put(newname, new TreeMap<String, String>());
                }

                for (String member : members) {
                    String[] names = member.split("_");

                    String name = names[0];
                    String id = names[1];
                    String contig = seqMap.get(name).get(id);

                    newSeqMap.get(newname).put(id, contig);
                }
            }
        } catch (FileNotFoundException e) {
            throw new IndianaException("Unable to cluster seqs", e);
        }

        return newSeqMap;
    }

    public void initialize() {
        log.info("Processing classifications...");
        maxContigs = CLASSIFICATIONS.size();

        List<ReferenceSequence> seqs = new ArrayList<ReferenceSequence>();

        for (String id : CLASSIFICATIONS.keySet()) {
            TableReader tr = new TableReader(CLASSIFICATIONS.get(id));

            int found = 0;

            for (Map<String, String> te : tr) {
                int kmersHB3 = Integer.valueOf(te.get("HB3"));
                int kmers3D7 = Integer.valueOf(te.get("3D7"));

                if (kmersHB3 > 1 && kmers3D7 > 1) {
                    String name = te.get("name");
                    String contig = te.get("contig");

                    if (id.equals("pe2")) {
                        if (!nameIdContigMap.containsKey(name)) {
                            nameIdContigMap.put(name, new TreeMap<String, String>());
                            nameIdPosMap.put(name, new HashMap<String, Integer>());
                        }
                        nameIdContigMap.get(name).put(id, contig);
                        nameIdPosMap.get(name).put(id, 0);
                    } else {
                        ReferenceSequence seq = new ReferenceSequence(id, 0, contig.getBytes());
                        seqs.add(seq);
                    }

                    found++;

                    if (contig.length() > longestContig) {
                        longestContig = contig.length();
                    }
                }
            }

            log.info("  Found {} chimeric contigs in {}", found, id);
        }

        log.info("  Found {} chimeric contigs total", nameIdContigMap.size() + seqs.size());

        for (String name : nameIdContigMap.keySet()) {
            String contig = nameIdContigMap.get(name).get("pe2");

            for (ReferenceSequence rseq : seqs) {
                String seq = new String(rseq.getBases());

                if (contig.contains(seq)) {
                    nameIdContigMap.get(name).put(rseq.getName(), seq);
                    nameIdPosMap.get(name).put(rseq.getName(), contig.indexOf(seq));
                }
            }
        }

        /*
        nameIdContigMap = clusterSequences(nameIdContigMap);
        log.info("  Grouped chimeric contigs into {} bins", nameIdContigMap.size());

        for (String name : nameIdContigMap.keySet()) {
            String contig = null;
            if      (nameIdContigMap.get(name).containsKey("pe2")) { contig = nameIdContigMap.get(name).get("pe2"); }
            else if (nameIdContigMap.get(name).containsKey("pe1")) { contig = nameIdContigMap.get(name).get("pe1"); }
            else if (nameIdContigMap.get(name).containsKey("se"))  { contig = nameIdContigMap.get(name).get("se"); }
            else if (nameIdContigMap.get(name).containsKey("sn"))  { contig = nameIdContigMap.get(name).get("sn"); }

            if (contig != null) {
                for (String id : nameIdContigMap.get(name).keySet()) {
                    if (!nameIdPosMap.containsKey(name)) {
                        nameIdPosMap.put(name, new HashMap<String, Integer>());
                        nameIdPosMap.get(name).put(id, 0);
                    }

                    String acontig = nameIdContigMap.get(name).get(id);

                    if (contig.contains(acontig)) {
                        nameIdPosMap.get(name).put(id, contig.indexOf(acontig));
                    }
                }

                log.info("{}: {}", name, nameIdContigMap.get(name));
            }
        }
        */

        //log.info("{}", nameIdContigMap);

        log.info("Processing references...");
        Color[] colors = generateColors(REFERENCES.size());

        int colorIndex = 0;
        for (String refid : REFERENCES.keySet()) {
            log.info("  {}:{}", refid, colors[colorIndex]);

            IndexedFastaSequenceFile ref = REFERENCES.get(refid);

            ReferenceSequence rseq;
            while ((rseq = ref.nextSequence()) != null) {
                String seq = new String(rseq.getBases());

                for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                    CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                    if (!kmerColorMap.containsKey(kmer)) {
                        kmerColorMap.put(kmer, colors[colorIndex]);
                    } else {
                        kmerColorMap.put(kmer, Color.GRAY);
                    }
                }
            }

            colorIndex++;
        }
    }

    private Color[] generateColors(int n) {
        Color[] cols = new Color[n];

        for(int i = 0; i < n; i++) {
            cols[i] = Color.getHSBColor((float) i / (float) n, 0.85f, 1.0f);
        }

        return cols;
    }

    public void setup() {
        size(marginLabel + 2*marginX + longestContig + marginXRight, marginTitle + 2*marginY + (marginContig + contigHeight)*maxContigs, PGraphicsPDF.PDF, out.getAbsolutePath());

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
        order.add("sn");
        order.add("se");

        log.info("Making plots...");

        int index = 0;
        for (String name : nameIdContigMap.keySet()) {
            log.info("  page {}", name);
            if (index > 0) { pdf.nextPage(); }
            index++;

            background(Color.WHITE.getRGB());

            fill(Color.BLACK.getRGB());
            textAlign(LEFT, CENTER);
            textSize(15);
            text(name, marginX + 30, marginTitle / 2);

            int i = 0;
            for (String id : order) {
                if (nameIdContigMap.get(name).containsKey(id)) {
                    String contig = nameIdContigMap.get(name).get(id);

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
                    text(contig.length() + " bp", xpos + contig.length() + 3, ypos - 2);

                    for (int k = 0; k <= contig.length() - KMER_SIZE; k++) {
                        CortexKmer kmer = new CortexKmer(contig.substring(k, k + KMER_SIZE));

                        Color color = kmerColorMap.containsKey(kmer) ? kmerColorMap.get(kmer) : Color.WHITE;

                        stroke(color.getRGB(), 150.0f);
                        strokeCap(PROJECT);
                        line(xpos + k, ypos + 1, xpos + k, ypos + contigHeight - 1);
                    }

                    stroke(Color.BLACK.getRGB());
                    noFill();
                    rect(xpos - 1, ypos, contig.length() + 1, contigHeight);

                    i++;
                }
            }
        }
    }

    public void draw() {
        exit();
    }
}
