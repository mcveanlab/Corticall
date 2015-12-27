package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import uk.ac.ox.well.indiana.Indiana;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.BwaAligner;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.LastzAligner;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataTables;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

public class GenotypeGraph extends Module {
    @Argument(fullName = "graphClean", shortName = "c", doc = "Graph (clean)")
    public CortexGraph CLEAN;

    @Argument(fullName = "graphDirty", shortName = "d", doc = "Graph (dirty)", required = false)
    public CortexGraph DIRTY;

    @Argument(fullName = "novelGraph", shortName = "n", doc = "Graph of novel kmers")
    public CortexGraph NOVEL;

    @Argument(fullName = "ref", shortName = "r", doc = "Fasta file for finished reference sequence")
    public File REF;

    @Argument(fullName = "ref1", shortName = "r1", doc = "Fasta file for first parent")
    public File REF1;

    @Argument(fullName = "ref2", shortName = "r2", doc = "Fasta file for second parent")
    public File REF2;

    @Argument(fullName = "bed", shortName = "b", doc = "Bed file describing variants", required = false)
    public File BED;

    @Argument(fullName = "novelKmerMap", shortName = "m", doc = "Novel kmer map", required = false)
    public File NOVEL_KMER_MAP;

    @Argument(fullName = "skipToKmer", shortName = "s", doc = "Skip processing to given kmer", required = false)
    public String KMER;

    @Argument(fullName = "haplotypes", shortName="hap", doc="Haplotypes file", required=false)
    public File HAPLOTYPES;

    @Argument(fullName = "novelKmerLimit", shortName="l", doc="Novel kmer count limit")
    public Integer NOVEL_KMER_LIMIT = 10000;

    @Output
    public PrintStream out;

    @Output(fullName="cout", shortName="co", doc="Circos track output")
    public PrintStream cout;

    private void evalVariant(GraphicalVariantContext gvc, int color, Map<CortexKmer, VariantInfo> vis, String stretch) {
        Set<VariantInfo> relevantVis = new HashSet<VariantInfo>();
        int kmerSize = vis.keySet().iterator().next().length();

        for (int i = 0; i <= stretch.length() - kmerSize; i++) {
            CortexKmer ck = new CortexKmer(stretch.substring(i, i + kmerSize));

            if (vis.containsKey(ck)) {
                relevantVis.add(vis.get(ck));
            }
        }

        if (relevantVis.size() > 0) {
            VariantInfo bestVi = null;

            for (VariantInfo vi : relevantVis) {
                String ref = gvc.getAttributeAsString(color, "parentalAllele");
                String alt = gvc.getAttributeAsString(color, "childAllele");

                String refStretch = gvc.getAttributeAsString(color, "parentalStretch");
                String altStretch = gvc.getAttributeAsString(color, "childStretch");

                int pos = gvc.getAttributeAsInt(color, "start");
                int refLength = gvc.getAttributeAsString(color, "parentalAllele").length();
                int altLength = gvc.getAttributeAsString(color, "childAllele").length();
                boolean found = false;

                while (pos >= 0 && pos + refLength < refStretch.length() && pos + altLength < altStretch.length()) {
                    String refFw = refStretch.substring(pos, pos + refLength);
                    String refRc = SequenceUtils.reverseComplement(refFw);

                    String altFw = altStretch.substring(pos, pos + altLength);
                    String altRc = SequenceUtils.reverseComplement(altFw);

                    if (vi.ref.equals(refFw) && vi.alt != null && vi.alt.equals(altFw)) {
                        ref = refFw;
                        alt = altFw;
                        found = true;
                        break;
                    } else if (vi.ref.equals(refRc) && vi.alt != null && vi.alt.equals(altRc)) {
                        ref = refRc;
                        alt = altRc;
                        found = true;
                        break;
                    }

                    pos--;
                }

                if (!found) {
                    pos = gvc.getAttributeAsInt(color, "start");

                    while (pos >= 0 && pos + refLength < refStretch.length() && pos + altLength < altStretch.length()) {
                        String refFw = refStretch.substring(pos, pos + refLength);
                        String refRc = SequenceUtils.reverseComplement(refFw);

                        String altFw = altStretch.substring(pos, pos + altLength);
                        String altRc = SequenceUtils.reverseComplement(altFw);

                        if (vi.ref.equals(refFw) && vi.alt != null && vi.alt.equals(altFw)) {
                            ref = refFw;
                            alt = altFw;
                            found = true;
                            break;
                        } else if (vi.ref.equals(refRc) && vi.alt != null && vi.alt.equals(altRc)) {
                            ref = refRc;
                            alt = altRc;
                            found = true;
                            break;
                        }

                        pos++;
                    }
                }

                String knownRef = vi.ref;
                String knownAlt = vi.alt == null ? "" : vi.alt;

                if (!found && knownRef.length() > ref.length() && knownAlt.length() > alt.length() && knownRef.length() == knownAlt.length()) {
                    String refFw = ref;
                    String altFw = alt;
                    String refRc = SequenceUtils.reverseComplement(refFw);
                    String altRc = SequenceUtils.reverseComplement(altFw);

                    String refFinal = null, altFinal = null;

                    if (knownRef.contains(refFw) && knownAlt.contains(altFw)) {
                        refFinal = refFw;
                        altFinal = altFw;
                    } else if (knownRef.contains(refRc) && knownAlt.contains(altRc)) {
                        refFinal = refRc;
                        altFinal = altRc;
                    }

                    if (refFinal != null && altFinal != null) {
                        int r0index = knownRef.indexOf(refFinal);
                        int a0index = knownAlt.indexOf(altFinal);
                        int r1index = r0index + refFinal.length();
                        int a1index = a0index + altFinal.length();

                        if (r0index == a0index && r1index == a1index &&
                                knownRef.substring(0, r0index).equals(knownAlt.substring(0, a0index)) &&
                                knownRef.substring(r1index, knownRef.length()).equals(knownAlt.substring(a1index, knownAlt.length()))) {
                            knownRef = ref;
                            knownAlt = alt;
                        }
                    }
                }

                log.debug("    - vi: {}", vi);
                log.debug("        known ref: {}", knownRef);
                log.debug("        known alt: {}", knownAlt);
                log.debug("        called ref: {}", ref);
                log.debug("        called alt: {}", alt);

                log.debug("    - matches: {} ({} {} {} {})", knownRef.equals(ref) && knownAlt.equals(alt), ref, alt, knownRef, knownAlt);

                vi.found = true;
                vi.matches = (knownRef.equals(ref) && knownAlt.equals(alt)) || vi.denovo.equals(gvc.getAttributeAsString(color, "event"));

                if (bestVi == null || vi.matches) {
                    bestVi = vi;
                }
            }

            if (bestVi != null) {
                String knownRef = bestVi.ref;
                String knownAlt = bestVi.alt != null ? bestVi.alt : "";
                String ref = gvc.getAttributeAsString(color, "parentalAllele");
                String alt = gvc.getAttributeAsString(color, "childAllele");

                gvc.attribute(color, "isKnownVariant", true);
                gvc.attribute(color, "variantId", bestVi.variantId);
                gvc.attribute(color, "allelesMatch", (knownRef.equals(ref) && knownAlt.equals(alt)) || bestVi.denovo.equals(gvc.getAttributeAsString(color, "event")));
                gvc.attribute(color, "eventsMatch", bestVi.denovo.equals(gvc.getAttributeAsString(color, "event")));
                gvc.attribute(color, "knownRef", knownRef);
                gvc.attribute(color, "knownAlt", knownAlt);
            } else {
                gvc.attribute(color, "isKnownVariant", false);
                gvc.attribute(color, "variantId", "unknown");
                gvc.attribute(color, "allelesMatch", false);
                gvc.attribute(color, "eventsMatch", false);
                gvc.attribute(color, "knownRef", "unknown");
                gvc.attribute(color, "knownAlt", "unknown");
            }
        }
    }

    private IntervalTreeMap<Integer> loadHaplotypes() {
        TableReader tr = new TableReader(HAPLOTYPES);

        String currentChr = "none";
        int currentParent = -1;
        int startPos = 0;
        int endPos = 0;

        IntervalTreeMap<Integer> itm = new IntervalTreeMap<Integer>();

        for (Map<String, String> te : tr) {
            for (String columnName : te.keySet()) {
                if (columnName.contains("/")) {
                    String[] pieces = columnName.split("/");
                    String newColumnName = pieces[1] + "." + pieces[2];

                    if (CLEAN.getColor(0).getSampleName().contains(newColumnName)) {
                        String chr = te.get("CHROM");
                        endPos = Integer.valueOf(te.get("POS"));
                        int parent = Integer.valueOf(te.get(columnName));

                        if (!currentChr.equals(chr)) {
                            startPos = 0;
                        }

                        if (currentParent != parent || !currentChr.equals(chr)) {
                            if (!currentChr.equals("none")) {
                                Interval interval = new Interval(currentChr, startPos, endPos);

                                itm.put(interval, currentParent);
                            }

                            currentChr = chr;
                            currentParent = parent;
                            startPos = endPos;
                        }
                    }
                }
            }
        }

        Interval interval = new Interval(currentChr, startPos, endPos);

        itm.put(interval, currentParent);

        return itm;
    }

    private int closestHaplotypeBackground(int[] hb, int index) {
        int hbFinal = 0;
        if (hb[index] == 1 || hb[index] == 2) {
            hbFinal = hb[index];
        } else if (hb[index] == 3) {
            int hbp = 3, dp = Integer.MAX_VALUE, hbn = 3, dn = Integer.MAX_VALUE;

            for (int i = index; i < hb.length; i++) {
                if (hb[i] == 1 || hb[i] == 2) {
                    hbn = hb[i];
                    dn = i - index;
                }
            }

            for (int i = index; i >= 0; i--) {
                if (hb[i] == 1 || hb[i] == 2) {
                    hbp = hb[i];
                    dp = index - i;
                }
            }

            hbFinal = (dp < dn) ? hbp : hbn;
        }

        return hbFinal;
    }

    private int[] getHaplotypeBackgroundVector(String alignmentStretch, KmerLookup kl, KmerLookup kl1, KmerLookup kl2, IntervalTreeMap<Integer> itm) {
        int[] hb = new int[alignmentStretch.length() - CLEAN.getKmerSize() + 1];

        for (int i = 0; i <= alignmentStretch.length() - CLEAN.getKmerSize(); i++) {
            String sk = alignmentStretch.substring(i, i + CLEAN.getKmerSize());
            CortexKmer ck = new CortexKmer(sk);
            CortexRecord cr = CLEAN.findRecord(ck);
            if (cr == null) {
                cr = DIRTY.findRecord(ck);
            }

            //log.info("  {} {} {}", sk, ck, cr);

            if (cr != null) {
                if      (cr.getCoverage(1) >  0 && cr.getCoverage(2) == 0) { hb[i] = 1; }
                else if (cr.getCoverage(1) == 0 && cr.getCoverage(2) >  0) { hb[i] = 2; }
                else if (cr.getCoverage(1) >  0 && cr.getCoverage(2) >  0) { hb[i] = 3; }
                else { hb[i] = 0; }
            }
        }

        for (int i = 0; i < hb.length; i++) {
            int oldhb = hb[i];

            if (hb[i] == 3) {
                hb[i] = closestHaplotypeBackground(hb, i);
            }
            int newhb = hb[i];

            if (hb[i] == 3) {
                String sk = alignmentStretch.substring(i, i + CLEAN.getKmerSize());

                Set<Interval> sis = kl.findKmer(sk);
                if (sis.size() == 0) { sis = kl.findKmer(SequenceUtils.reverseComplement(sk)); }

                if (sis.size() == 1) {
                    Interval si = sis.iterator().next();

                    if (itm.containsOverlapping(si)) {
                        hb[i] = itm.getOverlapping(si).iterator().next();
                    }
                } else {
                    sis = kl1.findKmer(sk);
                    if (sis.size() == 0) { sis = kl1.findKmer(SequenceUtils.reverseComplement(sk)); }

                    if (sis.size() == 1) {
                        hb[i] = 1;
                    } else {
                        sis = kl2.findKmer(sk);
                        if (sis.size() == 0) { sis = kl2.findKmer(SequenceUtils.reverseComplement(sk)); }

                        if (sis.size() == 1) {
                            hb[i] = 2;
                        }
                    }
                }
            }
        }

        for (int i = 0; i < hb.length; i++) {
            if (hb[i] == 3) {
                hb[i] = closestHaplotypeBackground(hb, i);
            }
        }

//        for (int i = 0; i < hb.length; i++) {
//            log.info("  {}", hb[i]);
//        }

        return hb;
    }

    private List<SAMRecord> getAlignment(String alignmentStretch, KmerLookup kl, KmerLookup kl1, KmerLookup kl2, IntervalTreeMap<Integer> itm) {
        int[] hb = getHaplotypeBackgroundVector(alignmentStretch, kl, kl1, kl2, itm);

        Map<String, Integer> pieces = new LinkedHashMap<String, Integer>();
        int currentHb = hb[0];
        int currentStart = 0;
        for (int i = 1; i < hb.length; i++) {
            if (hb[i] != currentHb) {
                String piece = alignmentStretch.substring(currentStart, i + CLEAN.getKmerSize());
                pieces.put(piece, currentHb);

                currentHb = hb[i];
                currentStart = i;
            }
        }

        String piece = alignmentStretch.substring(currentStart, alignmentStretch.length());
        pieces.put(piece, hb[currentStart]);

        Map<String, SAMRecord> alignments = new LinkedHashMap<String, SAMRecord>();
        Map<String, Integer> parentage = new LinkedHashMap<String, Integer>();

        LastzAligner la = new LastzAligner();
        for (String p : pieces.keySet()) {
            List<SAMRecord> srs;
            if (pieces.get(p) == 1) {
                srs = la.align(p, REF1);
            } else if (pieces.get(p) == 2) {
                srs = la.align(p, REF2);
            } else {
                srs = new ArrayList<SAMRecord>();
            }

            if (srs.size() == 1) {
                alignments.put(p, srs.iterator().next());
            } else {
                alignments.put(p, null);
            }

            parentage.put(p, pieces.get(p));
        }

        List<SAMRecord> alignmentsSimplified = new ArrayList<SAMRecord>();

        for (String p : alignments.keySet()) {
            SAMRecord sr = alignments.get(p);
            int pa = parentage.get(p);

            if (sr != null) {
                sr.setAttribute("PA", pa);

                alignmentsSimplified.add(sr);
            }
        }

        return alignmentsSimplified;
    }

    @Override
    public void execute() {
        log.info("Loading reference indices for fast kmer lookup...");
        KmerLookup kl  = new KmerLookup(REF);
        KmerLookup kl1 = new KmerLookup(REF1);
        KmerLookup kl2 = new KmerLookup(REF2);

        IntervalTreeMap<Integer> itm = (HAPLOTYPES != null && HAPLOTYPES.exists()) ? loadHaplotypes() : new IntervalTreeMap<Integer>();

        Map<CortexKmer, VariantInfo> vis = GenotypeGraphUtils.loadNovelKmerMap(NOVEL_KMER_MAP, BED);
        Map<String, VariantInfo> vids = new HashMap<String, VariantInfo>();
        Map<String, Boolean> viSeen = new HashMap<String, Boolean>();
        for (VariantInfo vi : vis.values()) {
            vids.put(vi.variantId, vi);
            viSeen.put(vi.variantId, false);
        }

        Map<CortexKmer, Boolean> novelKmers = new HashMap<CortexKmer, Boolean>();
        for (CortexRecord cr : NOVEL) {
            novelKmers.put(new CortexKmer(cr.getKmerAsString()), true);
        }

        if (novelKmers.size() > NOVEL_KMER_LIMIT) {
            throw new IndianaException("Too many novel kmers (for now)");
        }

        int totalNovelKmersUsed = 0;
        int stretchNum = 1;

        DataTables evalTables = new DataTables();

        evalTables.addTable("variantCalls", "Table containing variant information");

        evalTables.addTable("variantStats", "Statistics on variants", "knownVariantId", "knownVariantEvent", "knownVariantLength", "variantId", "variantEvent", "variantLength", "numVertices", "novelKmersContained", "novelKmersUsed", "filterStatus");

        evalTables.addTable("discoveryStats", "Statistics on variant discovery", "tp", "tn", "fp", "fn");
        evalTables.getTable("discoveryStats").set("dummy", "tp", 0l);
        evalTables.getTable("discoveryStats").set("dummy", "tn", 0l);
        evalTables.getTable("discoveryStats").set("dummy", "fp", 0l);
        evalTables.getTable("discoveryStats").set("dummy", "fn", 0l);

        evalTables.addTable("alleleMatchStats", "Statistics on allele matches", "match", "mismatch");
        evalTables.getTable("alleleMatchStats").set("dummy", "match", 0l);
        evalTables.getTable("alleleMatchStats").set("dummy", "mismatch", 0l);

        evalTables.addTable("eventMatchStats", "Statistics on event matches", "match", "mismatch");
        evalTables.getTable("eventMatchStats").set("dummy", "match", 0l);
        evalTables.getTable("eventMatchStats").set("dummy", "mismatch", 0l);

        Set<CortexKmer> novelKmersToVisit = new HashSet<CortexKmer>();

        if (KMER != null) {
            novelKmersToVisit.add(new CortexKmer(KMER));
        } else {
            novelKmersToVisit.addAll(novelKmers.keySet());
        }

        int variantsMissingKmers = 0;

        log.info("Genotyping novel kmer stretches in graph...");
        Set<GraphicalVariantContext> gvcs = new LinkedHashSet<GraphicalVariantContext>();
        for (CortexKmer novelKmer : novelKmersToVisit) {
            if (novelKmers.get(novelKmer)) {
                // Walk the graph left and right of novelKmer and extract a novel stretch
                String stretch;
                if (DIRTY == null) {
                    stretch = CortexUtils.getSeededStretch(CLEAN, novelKmer.getKmerAsString(), 0, true);
                } else {
                    stretch = CortexUtils.getSeededStretch(CLEAN, DIRTY, novelKmer.getKmerAsString(), 0, true);
                }

                log.info("  stretch {}: {} bp", stretchNum, stretch.length());
                log.info("    novel kmer: {}", novelKmer);
                log.info("    sequence:");
                log.info("    - 0: {}", SequenceUtils.truncate(stretch, 100));

                // Construct GVC
                GraphicalVariantContext gvc = new GraphicalVariantContext()
                        .attribute(0, "stretch", stretch)
                        .attribute(0, "stretchNum", stretchNum)
                        .attribute(0, "stretchLength", stretch.length())
                        .attribute(0, "novelKmersTotal", novelKmers.size());

                // Fetch the local subgraph context from disk
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> ag = GenotypeGraphUtils.loadLocalSubgraph(stretch, CLEAN, DIRTY, novelKmers);
                int numPredecessors = 0, numSuccessors = 0;
                for (AnnotatedVertex av : ag.vertexSet()) {
                    if (av.flagIsSet("predecessor")) { numPredecessors++; }
                    if (av.flagIsSet("successor")) { numSuccessors++; }
                }

                log.info("    subgraph : {} vertices, {} edges, {} predecessors, {} successors", ag.vertexSet().size(), ag.edgeSet().size(), numPredecessors, numSuccessors);
                log.info("    alignment info...");

                GenotypeGraphUtils.annotateAlignmentInformation(ag, kl1, kl2);

                // Extract parental stretches
                log.info("    min weight paths...");
                PathInfo p1 = GenotypeGraphUtils.computeBestMinWeightPath(CLEAN, DIRTY, ag, 1, stretch, novelKmers);
                PathInfo p2 = GenotypeGraphUtils.computeBestMinWeightPath(CLEAN, DIRTY, ag, 2, stretch, novelKmers);

                log.info("    paths:");
                log.info("    - 1: {} ({} bp)", SequenceUtils.truncate(p1.parent, 100), p1.parent.length());
                log.info("      c: {} ({} bp)", SequenceUtils.truncate(p1.child, 100), p1.child.length());
                log.info("    - 2: {} ({} bp)", SequenceUtils.truncate(p2.parent, 100), p2.parent.length());
                log.info("      c: {} ({} bp)", SequenceUtils.truncate(p2.child, 100), p2.child.length());

                // Call variants
                log.info("    variants...");
                gvc.add(GenotypeGraphUtils.callVariant(CLEAN, DIRTY, p1, 1, stretch, ag, kl1));
                gvc.add(GenotypeGraphUtils.callVariant(CLEAN, DIRTY, p2, 2, stretch, ag, kl2));

                log.info("    variants:");
                log.info("    - 1: {} {} ({} bp)", gvc.getAttributeAsString(1, "event"), SequenceUtils.truncate(gvc.getAttributeAsString(1, "parentalAllele"), 70), gvc.getAttributeAsString(1, "parentalAllele").length());
                log.info("      c: {} {} ({} bp)", gvc.getAttributeAsString(1, "event"), SequenceUtils.truncate(gvc.getAttributeAsString(1, "childAllele"), 70), gvc.getAttributeAsString(1, "childAllele").length());
                log.info("    - 2: {} {} ({} bp)", gvc.getAttributeAsString(2, "event"), SequenceUtils.truncate(gvc.getAttributeAsString(2, "parentalAllele"), 70), gvc.getAttributeAsString(2, "parentalAllele").length());
                log.info("      c: {} {} ({} bp)", gvc.getAttributeAsString(2, "event"), SequenceUtils.truncate(gvc.getAttributeAsString(2, "childAllele"), 70), gvc.getAttributeAsString(2, "childAllele").length());

                // Finalize into a single call
                GenotypeGraphUtils.chooseVariant(gvc);

                // Show alignment
                log.info("    alignment:");
                log.info("    - novel stretch: {}", gvc.getAttribute(0, "novelStretchAlignment"));
                log.info("    - parental path: {}", gvc.getAttribute(0, "parentalPathAlignment"));

                // See how many novel kmers we've used up
                int novelKmersContained = 0;
                int novelKmersUsed = 0;

                for (AnnotatedVertex av : ag.vertexSet()) {
                    CortexKmer ck = new CortexKmer(av.getKmer());

                    if (novelKmers.containsKey(ck)) {
                        novelKmersContained++;

                        if (novelKmers.get(ck)) {
                            totalNovelKmersUsed++;
                            novelKmersUsed++;
                            novelKmers.put(ck, false);
                        }
                    }
                }

                gvc.attribute(0, "novelKmersContained", novelKmersContained);
                gvc.attribute(0, "novelKmersUsed", novelKmersUsed);

                log.info("    novelty:");
                log.info("    - novel kmers present:  {}/{}", novelKmersContained, novelKmers.size());
                log.info("    - novel kmers utilized: {}/{}", novelKmersUsed, novelKmers.size());
                log.info("    - cumulative usage:     {}/{}", totalNovelKmersUsed, novelKmers.size());

                gvc.attribute(0, "filter", (novelKmersContained <= 1 && gvc.getAttributeAsString(0, "event").equals("unknown")) ? "FAIL" : "PASS");
                log.info("    filter: {}", gvc.getAttributeAsString(0, "filter"));

                // Evaluate variants
                if (BED != null) {
                    log.info("    evaluate:");

                    for (int c = 0; c < 3; c++) {
                        evalVariant(gvc, c, vis, stretch);

                        String vid = gvc.getAttributeAsString(c, "variantId");
                        VariantInfo vi = vids.containsKey(vid) ? vids.get(vid) : null;
                        String event = vi == null ? "none" : vi.denovo;

                        log.info("    - {}: background {}, isKnownVariant {}, allelesMatch {}, eventsMatch {}, variantId {} {} {} {}",
                                c,
                                gvc.getAttributeAsInt(c, "haplotypeBackground"),
                                gvc.getAttribute(c, "isKnownVariant"),
                                gvc.getAttribute(c, "allelesMatch"),
                                gvc.getAttribute(c, "eventsMatch"),
                                gvc.getAttributeAsString(c, "variantId"),
                                event,
                                gvc.getAttributeAsString(c, "knownRef"),
                                gvc.getAttributeAsString(c, "knownAlt")
                        );
                    }

                    if (gvc.getAttributeAsBoolean(0, "isKnownVariant")) {
                        evalTables.getTable("discoveryStats").increment("dummy", "tp");

                        String vid = gvc.getAttributeAsString(0, "variantId");
                        VariantInfo vi = vids.get(vid);

                        int refLength = vi.ref != null ? vi.ref.length() : 0;
                        int altLength = vi.alt != null ? vi.alt.length() : 0;

                        String pk = vi.variantId + "." + gvc.getAttributeAsInt(0, "stretchNum");

                        evalTables.getTable("variantStats").set(pk, "knownVariantId", gvc.getAttributeAsString(0, "variantId"));
                        evalTables.getTable("variantStats").set(pk, "knownVariantEvent", vi.denovo);
                        evalTables.getTable("variantStats").set(pk, "knownVariantLength", Math.abs(refLength - altLength));
                        evalTables.getTable("variantStats").set(pk, "variantId", gvc.getAttributeAsInt(0, "stretchNum"));
                        evalTables.getTable("variantStats").set(pk, "variantEvent", gvc.getAttributeAsString(0, "event"));
                        evalTables.getTable("variantStats").set(pk, "variantLength", Math.abs(gvc.getAttributeAsString(0, "parentalAllele").length() - gvc.getAttributeAsString(0, "childAllele").length()));
                        evalTables.getTable("variantStats").set(pk, "numVertices", ag.vertexSet().size());
                        evalTables.getTable("variantStats").set(pk, "novelKmersContained", novelKmersContained);
                        evalTables.getTable("variantStats").set(pk, "novelKmersUsed", novelKmersUsed);
                        evalTables.getTable("variantStats").set(pk, "filterStatus", gvc.getAttributeAsString(0, "filter"));
                        evalTables.getTable("variantStats").set(pk, "seedKmer", novelKmer.getKmerAsString());
                        evalTables.getTable("variantStats").set(pk, "numPredecessors", numPredecessors);
                        evalTables.getTable("variantStats").set(pk, "numSuccessors", numSuccessors);

                        viSeen.put(vi.variantId, true);

                        if (gvc.getAttributeAsBoolean(0, "allelesMatch")) {
                            evalTables.getTable("alleleMatchStats").increment("dummy", "match");
                        } else {
                            evalTables.getTable("alleleMatchStats").increment("dummy", "mismatch");
                        }

                        if (gvc.getAttributeAsBoolean(0, "eventsMatch")) {
                            evalTables.getTable("eventMatchStats").increment("dummy", "match");
                        } else {
                            evalTables.getTable("eventMatchStats").increment("dummy", "mismatch");
                        }
                    } else {
                        if (gvc.getAttributeAsString(0, "filter").equals("PASS")) {
                            evalTables.getTable("discoveryStats").increment("dummy", "fp");
                        } else {
                            evalTables.getTable("discoveryStats").increment("dummy", "tn");
                        }

                        String pk = "none." + gvc.getAttributeAsInt(0, "stretchNum");

                        evalTables.getTable("variantStats").set(pk, "knownVariantId", "none");
                        evalTables.getTable("variantStats").set(pk, "knownVariantEvent", "none");
                        evalTables.getTable("variantStats").set(pk, "knownVariantLength", 0);
                        evalTables.getTable("variantStats").set(pk, "variantId", gvc.getAttributeAsInt(0, "stretchNum"));
                        evalTables.getTable("variantStats").set(pk, "variantEvent", gvc.getAttributeAsString(0, "event"));
                        evalTables.getTable("variantStats").set(pk, "variantLength", Math.abs(gvc.getAttributeAsString(0, "parentalAllele").length() - gvc.getAttributeAsString(0, "childAllele").length()));
                        evalTables.getTable("variantStats").set(pk, "numVertices", ag.vertexSet().size());
                        evalTables.getTable("variantStats").set(pk, "novelKmersContained", novelKmersContained);
                        evalTables.getTable("variantStats").set(pk, "novelKmersUsed", novelKmersUsed);
                        evalTables.getTable("variantStats").set(pk, "filterStatus", gvc.getAttributeAsString(0, "filter"));
                        evalTables.getTable("variantStats").set(pk, "seedKmer", novelKmer.getKmerAsString());
                        evalTables.getTable("variantStats").set(pk, "numPredecessors", numPredecessors);
                        evalTables.getTable("variantStats").set(pk, "numSuccessors", numSuccessors);
                    }
                }

                log.info("");

                gvcs.add(gvc);

                evalTables.getTable("variantCalls").set(novelKmer.getKmerAsString(), "novelKmer", novelKmer);
                evalTables.getTable("variantCalls").set(novelKmer.getKmerAsString(), "stretchNum", stretchNum);
                evalTables.getTable("variantCalls").set(novelKmer.getKmerAsString(), "event", gvc.getAttributeAsString(0, "event"));
                evalTables.getTable("variantCalls").set(novelKmer.getKmerAsString(), "parentalAlleleLength", gvc.getAttributeAsString(0, "parentalAllele").length());
                evalTables.getTable("variantCalls").set(novelKmer.getKmerAsString(), "childAlleleLength", gvc.getAttributeAsString(0, "childAllele").length());
                evalTables.getTable("variantCalls").set(novelKmer.getKmerAsString(), "parentalAllele", gvc.getAttributeAsString(0, "parentalAllele"));
                evalTables.getTable("variantCalls").set(novelKmer.getKmerAsString(), "childAllele", gvc.getAttributeAsString(0, "childAllele"));
                evalTables.getTable("variantCalls").set(novelKmer.getKmerAsString(), "score", gvc.getAttributeAsInt(0, "score"));
                evalTables.getTable("variantCalls").set(novelKmer.getKmerAsString(), "haplotypicBackground", gvc.getAttributeAsInt(0, "haplotypicBackground"));
                evalTables.getTable("variantCalls").set(novelKmer.getKmerAsString(), "start", gvc.getAttributeAsInt(0, "start"));
                evalTables.getTable("variantCalls").set(novelKmer.getKmerAsString(), "stop", gvc.getAttributeAsInt(0, "stop"));
                evalTables.getTable("variantCalls").set(novelKmer.getKmerAsString(), "novelKmersContained", novelKmersContained);
                evalTables.getTable("variantCalls").set(novelKmer.getKmerAsString(), "novelKmersUsed", novelKmersUsed);
                evalTables.getTable("variantCalls").set(novelKmer.getKmerAsString(), "novelKmersTotal", novelKmers.size());
                evalTables.getTable("variantCalls").set(novelKmer.getKmerAsString(), "traversalStatus", gvc.getAttributeAsString(0, "traversalStatus"));
                evalTables.getTable("variantCalls").set(novelKmer.getKmerAsString(), "stretch", stretch);
                evalTables.getTable("variantCalls").set(novelKmer.getKmerAsString(), "parentalStretch", gvc.getAttributeAsString(0, "parentalStretch"));
                evalTables.getTable("variantCalls").set(novelKmer.getKmerAsString(), "childStretch", gvc.getAttributeAsString(0, "childStretch"));

                String alignmentStretch = (!gvc.getAttributeAsString(0, "parentalStretch").isEmpty()) ? gvc.getAttributeAsString(0, "parentalStretch") : gvc.getAttributeAsString(0, "stretch");

                List<SAMRecord> alignments = getAlignment(alignmentStretch, kl, kl1, kl2, itm);

                for (SAMRecord sr : alignments) {
                    log.info("  {}", sr.getSAMString());
                }

                if (alignments.size() == 1) {
                    SAMRecord sr = alignments.get(0);
                    int radius = sr.getIntegerAttribute("PA") == 1 ? 8 : 6;

                    if (gvc.getAttribute(0, "event").equals("GC")) {
                        cout.printf("stretch%d %s %d %d radius1=0.%dr # %s_%s\n", stretchNum, sr.getReferenceName(), sr.getAlignmentStart(), sr.getAlignmentEnd(), 8, gvc.getAttribute(0, "event"), gvc.getAttribute(0, "traversalStatus"));
                        cout.printf("stretch%d %s %d %d radius2=0.%dr # %s_%s\n", stretchNum, sr.getReferenceName(), sr.getAlignmentStart(), sr.getAlignmentEnd(), 6, gvc.getAttribute(0, "event"), gvc.getAttribute(0, "traversalStatus"));
                    } else {
                        cout.printf("stretch%d %s %d %d radius1=0.%dr # %s_%s\n", stretchNum, sr.getReferenceName(), sr.getAlignmentStart(), sr.getAlignmentEnd(), radius, gvc.getAttribute(0, "event"), gvc.getAttribute(0, "traversalStatus"));
                        cout.printf("stretch%d %s %d %d radius2=0.%dr # %s_%s\n", stretchNum, sr.getReferenceName(), sr.getAlignmentStart(), sr.getAlignmentEnd(), radius, gvc.getAttribute(0, "event"), gvc.getAttribute(0, "traversalStatus"));
                    }
                } else if (alignments.size() > 1) {
                    for (int i = 0; i < alignments.size() - 1; i++) {
                        SAMRecord sr1 = alignments.get(i);
                        SAMRecord sr2 = alignments.get(i + 1);

                        int r1 = sr1.getIntegerAttribute("PA") == 1 ? 8 : 6;
                        int r2 = sr2.getIntegerAttribute("PA") == 1 ? 8 : 6;

                        cout.printf("stretch%d.%d %s %d %d radius1=0.%dr # %s_%s\n", stretchNum, i, sr1.getReferenceName(), sr1.getAlignmentStart(), sr1.getAlignmentEnd(), r1, gvc.getAttribute(0, "event"), gvc.getAttribute(0, "traversalStatus"));
                        cout.printf("stretch%d.%d %s %d %d radius2=0.%dr # %s_%s\n", stretchNum, i, sr2.getReferenceName(), sr2.getAlignmentStart(), sr2.getAlignmentEnd(), r2, gvc.getAttribute(0, "event"), gvc.getAttribute(0, "traversalStatus"));
                    }
                } else {
                    cout.printf("stretch%d %s %d %d radius1=0.%dr # %s_%s\n", stretchNum, "NA", 100*stretchNum, 100*stretchNum, 7, gvc.getAttribute(0, "event"), gvc.getAttribute(0, "traversalStatus"));
                    cout.printf("stretch%d %s %d %d radius2=0.%dr # %s_%s\n", stretchNum, "NA", 100*stretchNum, 100*stretchNum, 7, gvc.getAttribute(0, "event"), gvc.getAttribute(0, "traversalStatus"));
                }

                stretchNum++;
            }
        }

        log.info("Variants missing kmers: {}", variantsMissingKmers);

        for (String vid : viSeen.keySet()) {
            if (!viSeen.get(vid)) {
                evalTables.getTable("discoveryStats").increment("dummy", "fn");

                evalTables.getTable("variantStats").set(vid, "knownVariantId", vid);
                evalTables.getTable("variantStats").set(vid, "knownVariantEvent", vids.get(vid).denovo);

                String pk = vid + ".none";

                VariantInfo vi = vids.get(vid);

                int refLength = vi.ref != null ? vi.ref.length() : 0;
                int altLength = vi.alt != null ? vi.alt.length() : 0;

                evalTables.getTable("variantStats").set(pk, "knownVariantId", vid);
                evalTables.getTable("variantStats").set(pk, "knownVariantEvent", vids.get(vid).denovo);
                evalTables.getTable("variantStats").set(pk, "knownVariantLength", Math.abs(refLength - altLength));
                evalTables.getTable("variantStats").set(pk, "variantId", "none");
                evalTables.getTable("variantStats").set(pk, "variantEvent", "none");
                evalTables.getTable("variantStats").set(pk, "variantLength", 0);
                evalTables.getTable("variantStats").set(pk, "numVertices", 0);
                evalTables.getTable("variantStats").set(pk, "novelKmersContained", 0);
                evalTables.getTable("variantStats").set(pk, "novelKmersUsed", 0);
                evalTables.getTable("variantStats").set(pk, "filterStatus", "NONE");
                evalTables.getTable("variantStats").set(pk, "seedKmer", "none");
                evalTables.getTable("variantStats").set(pk, "numPredecessors", 0);
                evalTables.getTable("variantStats").set(pk, "numSuccessors", 0);
            }
        }

        evalTables.write(out);
    }
}
