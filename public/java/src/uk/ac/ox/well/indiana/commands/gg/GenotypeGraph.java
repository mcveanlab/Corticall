package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import uk.ac.ox.well.indiana.Indiana;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.BwaAligner;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.ExternalAligner;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.LastzAligner;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.LocalAligner;
import uk.ac.ox.well.indiana.utils.alignment.sw.SmithWaterman;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataTable;
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

    @Argument(fullName = "rejectGraph", shortName = "x", doc = "Graph of reject kmers")
    public CortexGraph REJECT;

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
                    String newColumnName = pieces.length > 2 ? pieces[1] + "." + pieces[2] : pieces[0] + "." + pieces[1];

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

        ExternalAligner la = new LastzAligner();
        //ExternalAligner la = new BwaAligner();
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

    private void smooth(List<List<Interval>> al, String tag) {
        DataTable dt = new DataTable(tag, "alignment");

        for (int i = 0; i < al.size(); i++) {
            if (al.get(i).size() == 0) {
                dt.set("0", String.valueOf(i), "none");
            } else {
                for (int j = 0; j < al.get(i).size(); j++) {
                    dt.set(String.valueOf(j), String.valueOf(i), String.format("%s:%d", al.get(i).get(j).getSequence(), al.get(i).get(j).getStart()));
                }
            }
        }

        log.info("\n{}", dt);
    }

    private StringBuilder smooth(StringBuilder sb) {
        for (int i = 0; i < sb.length(); i++) {
            if (sb.charAt(i) == 'B' || sb.charAt(i) == '?') {
                char leftChar = '?', rightChar = '?';
                int leftDist = Integer.MAX_VALUE, rightDist = Integer.MAX_VALUE;

                for (int j = i - 1; j >= 0; j--) {
                    if (sb.charAt(j) == '1' || sb.charAt(j) == '2') {
                        leftChar = sb.charAt(j);
                        leftDist = i - j;
                    }
                }

                for (int j = i + 1; j < sb.length(); j++) {
                    if (sb.charAt(j) == '1' || sb.charAt(j) == '2') {
                        rightChar = sb.charAt(j);
                        rightDist = j - i;
                    }
                }

                char bestChr = leftDist <= rightDist ? leftChar : rightChar;
                sb.setCharAt(i, bestChr);
            }
        }

        return sb;
    }

    private List<Interval> combineIntervals(List<List<Interval>> allIntervals) {
        List<Interval> combinedIntervals = new ArrayList<Interval>();

        Interval currentInterval = null;
        for (int i = 0; i < allIntervals.size(); i++) {
            List<Interval> intervals = allIntervals.get(i);
            Interval interval;

            if (intervals.size() == 1) {
                interval = intervals.iterator().next();
            } else {
                interval = new Interval("*", i, i + 1);
            }

            if (currentInterval == null) {
                currentInterval = new Interval(interval.getSequence(), interval.getStart(), interval.getEnd(), interval.isNegativeStrand(), "none");
            } else {
                if (currentInterval.getSequence().equals(interval.getSequence()) && interval.intersects(currentInterval)) {
                    currentInterval = new Interval(
                            currentInterval.getSequence(),
                            currentInterval.getStart() < interval.getStart() ? currentInterval.getStart() : interval.getStart(),
                            currentInterval.getEnd()   > interval.getEnd()   ? currentInterval.getEnd()   : interval.getEnd(),
                            interval.isNegativeStrand(),
                            "."
                    );
                } else {
                    combinedIntervals.add(currentInterval);

                    currentInterval = interval;
                }
            }
        }

        if (currentInterval != null) {
            combinedIntervals.add(currentInterval);
        }

        return combinedIntervals;
    }

    private Set<CortexKmer> loadRejectedKmers() {
        Set<CortexKmer> rejects = new HashSet<CortexKmer>();

        for (CortexRecord cr : REJECT) {
            rejects.add(cr.getCortexKmer());
        }

        return rejects;
    }

    private String compactAlignments(List<Interval> intervals) {
        List<String> compactIntervals = new ArrayList<String>();

        for (Interval interval : intervals) {
            compactIntervals.add(String.format("%s:%d-%d:%s", interval.getSequence(), interval.getStart(), interval.getEnd(), interval.isPositiveStrand() ? "+" : "-"));
        }

        return Joiner.on(";").join(compactIntervals);
    }

    @Override
    public void execute() {
        Random rnd = new Random(0);

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

        Set<CortexKmer> rejects = loadRejectedKmers();

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
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> ag = GenotypeGraphUtils.loadLocalSubgraph(stretch, CLEAN, DIRTY, novelKmers, false);
                int numPredecessors = 0, numSuccessors = 0;
                for (AnnotatedVertex av : ag.vertexSet()) {
                    if (av.flagIsSet("predecessor")) { numPredecessors++; }
                    if (av.flagIsSet("successor")) { numSuccessors++; }
                }

                log.info("    subgraph : {} vertices, {} edges, {} predecessors, {} successors", ag.vertexSet().size(), ag.edgeSet().size(), numPredecessors, numSuccessors);

                GenotypeGraphUtils.annotateAlignmentInformation(ag, kl1, kl2);

                // Extract parental stretches
                PathInfo p1 = GenotypeGraphUtils.computeBestMinWeightPath(CLEAN, DIRTY, ag, 1, stretch, novelKmers);
                PathInfo p2 = GenotypeGraphUtils.computeBestMinWeightPath(CLEAN, DIRTY, ag, 2, stretch, novelKmers);

                log.info("    - 1: {} ({} bp)", SequenceUtils.truncate(p1.parent, 100), p1.parent.length());
                log.info("      c: {} ({} bp)", SequenceUtils.truncate(p1.child, 100), p1.child.length());
                log.info("    - 2: {} ({} bp)", SequenceUtils.truncate(p2.parent, 100), p2.parent.length());
                log.info("      c: {} ({} bp)", SequenceUtils.truncate(p2.child, 100), p2.child.length());

                // Call variants
                gvc.add(GenotypeGraphUtils.callVariant(CLEAN, DIRTY, p1, 1, stretch, ag, kl1));
                gvc.add(GenotypeGraphUtils.callVariant(CLEAN, DIRTY, p2, 2, stretch, ag, kl2));

                log.info("    candidates:");
                log.info("    - 1: {} {} ({} bp)", gvc.getAttributeAsString(1, "event"), SequenceUtils.truncate(gvc.getAttributeAsString(1, "parentalAllele"), 80), gvc.getAttributeAsString(1, "parentalAllele").length());
                log.info("      c: {} {} ({} bp)", gvc.getAttributeAsString(1, "event"), SequenceUtils.truncate(gvc.getAttributeAsString(1, "childAllele"), 80), gvc.getAttributeAsString(1, "childAllele").length());
                log.info("    - 2: {} {} ({} bp)", gvc.getAttributeAsString(2, "event"), SequenceUtils.truncate(gvc.getAttributeAsString(2, "parentalAllele"), 80), gvc.getAttributeAsString(2, "parentalAllele").length());
                log.info("      c: {} {} ({} bp)", gvc.getAttributeAsString(2, "event"), SequenceUtils.truncate(gvc.getAttributeAsString(2, "childAllele"), 80), gvc.getAttributeAsString(2, "childAllele").length());

                // Finalize into a single call
                GenotypeGraphUtils.chooseVariant(gvc);

                int h = gvc.getAttributeAsInt(0, "haplotypeBackground") <= 1 ? 1 : 2;

                log.info("    variant:");
                log.info("    - {}: {} {} ({} bp)", h, gvc.getAttributeAsString(h, "event"), SequenceUtils.truncate(gvc.getAttributeAsString(h, "parentalAllele"), 80), gvc.getAttributeAsString(h, "parentalAllele").length());
                log.info("      c: {} {} ({} bp)",     gvc.getAttributeAsString(h, "event"), SequenceUtils.truncate(gvc.getAttributeAsString(h, "childAllele"), 80), gvc.getAttributeAsString(h, "childAllele").length());

                String pstretch = gvc.getAttributeAsString(0, "parentalStretch");
                String cstretch = gvc.getAttributeAsString(0, "childStretch");

                String astretch = (pstretch.isEmpty() && cstretch.isEmpty()) ? stretch : pstretch;
                String bstretch = (pstretch.isEmpty() && cstretch.isEmpty()) ? stretch : cstretch;

                List<Interval> finalPos = new ArrayList<Interval>();
                List<List<Interval>> alignment = gvc.getAttributeAsInt(0, "haplotypeBackground") == 1 ? kl1.alignSmoothly(astretch) : kl2.alignSmoothly(astretch);
                for (int i = 0; i < alignment.size(); i++) {
                    for (int j = 0; j < alignment.get(i).size(); j++) {
                        Interval oldInterval = alignment.get(i).get(j);
                        Interval newInterval = new Interval(oldInterval.getSequence(), oldInterval.getStart(), oldInterval.getEnd(), oldInterval.isNegativeStrand(), String.valueOf(gvc.getAttributeAsInt(0, "haplotypeBackground")));

                        alignment.get(i).set(j, newInterval);
                    }
                }

                if (pstretch.isEmpty() && cstretch.isEmpty()) {
                    finalPos = combineIntervals(alignment);

                    StringBuilder nv = new StringBuilder();
                    for (int i = 0; i <= astretch.length() - CLEAN.getKmerSize(); i++) {
                        CortexKmer ck = new CortexKmer(astretch.substring(i, i + CLEAN.getKmerSize()));
                        CortexRecord cr = CLEAN.findRecord(ck);
                        if (cr == null) {
                            cr = DIRTY.findRecord(ck);
                        }

                        if (cr != null) {
                            if (cr.getCoverage(0) > 0 && cr.getCoverage(1) == 0 && cr.getCoverage(2) == 0) {
                                nv.append(".");
                            } else if (cr.getCoverage(1) > 0 && cr.getCoverage(2) == 0) {
                                nv.append("1");
                            } else if (cr.getCoverage(1) == 0 && cr.getCoverage(2) > 0) {
                                nv.append("2");
                            } else if (cr.getCoverage(1) > 0 && cr.getCoverage(2) > 0) {
                                nv.append("B");
                            } else {
                                nv.append("?");
                            }
                        } else {
                            nv.append(" ");
                        }
                    }

                    nv = smooth(nv);

                    List<Pair<String, String>> sections = new ArrayList<Pair<String, String>>();
                    StringBuilder first = new StringBuilder();
                    StringBuilder second = new StringBuilder();
                    for (int i = 0; i < nv.length(); i++) {
                        if (nv.charAt(i) == '.' && (first.length() > 0 || second.length() > 0)) {
                            sections.add(new Pair<String, String>(first.toString(), second.toString()));
                            first = new StringBuilder();
                            second = new StringBuilder();
                        }

                        if (nv.charAt(i) == '1' || nv.charAt(i) == 'B') {
                            if (first.length() == 0) {
                                first.append(astretch.substring(i, i + CLEAN.getKmerSize()));
                            } else {
                                first.append(astretch.charAt(i + CLEAN.getKmerSize() - 1));
                            }
                        }

                        if (nv.charAt(i) == '2' || nv.charAt(i) == 'B') {
                            if (second.length() == 0) {
                                second.append(astretch.substring(i, i + CLEAN.getKmerSize()));
                            } else {
                                second.append(astretch.charAt(i + CLEAN.getKmerSize() - 1));
                            }
                        }
                    }

                    if (first.length() > 0 || second.length() > 0) {
                        sections.add(new Pair<String, String>(first.toString(), second.toString()));
                    }

                    ExternalAligner la = new LastzAligner();
                    for (Pair<String, String> p : sections) {
                        List<SAMRecord> afl = p.getFirst().isEmpty()  ? null : la.align(p.getFirst(),  REF1);
                        List<SAMRecord> asl = p.getSecond().isEmpty() ? null : la.align(p.getSecond(), REF2);

                        SAMRecord af = afl != null && afl.size() == 1 ? afl.get(0) : null;
                        SAMRecord as = asl != null && asl.size() == 1 ? asl.get(0) : null;
                        SAMRecord a = af;
                        int b = 0;

                        if (af != null && as == null) {
                            a = af;
                            b = 1;
                        } else if (af == null && as != null) {
                            a = as;
                            b = 2;
                        } else if (af != null && as != null) {
                            if (rnd.nextBoolean()) { a = af; b = 1; }
                            else { a = as; b = 2; }
                        }

                        if (a != null) {
                            finalPos.add(new Interval(a.getReferenceName(), a.getAlignmentStart(), a.getAlignmentEnd(), false, String.valueOf(b)));
                        }
                    }
                } else {
                    //smooth(alignment, "kl" + gvc.getAttributeAsInt(0, "haplotypeBackground"));

                    int startIndex = gvc.getAttributeAsInt(0, "start") - CLEAN.getKmerSize();
                    int endIndex = gvc.getAttributeAsInt(0, "e1") - CLEAN.getKmerSize();

                    Interval start = alignment.size() > startIndex && alignment.get(startIndex).size() > 0 ? alignment.get(startIndex).get(0) : null;
                    Interval end = alignment.size() > endIndex && alignment.get(endIndex).size() > 0 ? alignment.get(endIndex).get(0) : null;

                    if (start != null && end != null) {
                        Interval pos = new Interval(start.getSequence(), start.getStart() + startIndex + 1, start.getStart() + startIndex + 1);

                        if (start.getSequence().equals(end.getSequence())) {
                            if (start.isNegativeStrand() && end.isNegativeStrand()) {
                                pos = new Interval(end.getSequence(), end.getStart() + 1, end.getStart() + 1, true, String.valueOf(gvc.getAttributeAsInt(0, "haplotypeBackground")));
                            } else {
                                pos = new Interval(start.getSequence(), start.getStart() + 1, start.getStart() + 1, true, String.valueOf(gvc.getAttributeAsInt(0, "haplotypeBackground")));
                            }
                        }

                        finalPos.add(pos);

                        log.info("    {} {}", start, end);
                    } else if ((start == null && end != null) || (start != null && end == null)) {
                        Interval pos = null;

                        if (start != null) {
                            if (start.isNegativeStrand()) {
                                pos = new Interval(start.getSequence(), start.getStart() + 1, start.getStart() + 1, true, String.valueOf(gvc.getAttributeAsInt(0, "haplotypeBackground")));
                            } else {
                                pos = new Interval(start.getSequence(), start.getEnd() + 1, start.getEnd() + 1, false, String.valueOf(gvc.getAttributeAsInt(0, "haplotypeBackground")));
                            }
                        } else if (end != null) {
                            if (end.isNegativeStrand()) {
                                pos = new Interval(end.getSequence(), end.getStart() + 1, end.getStart() + 1, true, String.valueOf(gvc.getAttributeAsInt(0, "haplotypeBackground")));
                            } else {
                                pos = new Interval(end.getSequence(), end.getEnd() + 1, end.getEnd() + 1, false, String.valueOf(gvc.getAttributeAsInt(0, "haplotypeBackground")));
                            }
                        }

                        finalPos.add(pos);

                        log.info("    {} {}", start, end);
                    } else {
                        for (List<Interval> its : alignment) {
                            if (its.size() > 0) {
                                finalPos.add(its.get(0));
                                log.info("    {}", its.get(0));
                                break;
                            }
                        }
                    }
                }

                boolean hasDirtyKmers = false;
                for (int i = 0; i <= bstretch.length() - CLEAN.getKmerSize(); i++) {
                    CortexKmer ck = new CortexKmer(bstretch.substring(i, i + CLEAN.getKmerSize()));

                    if (CLEAN.findRecord(ck) == null && DIRTY.findRecord(ck) != null) {
                        hasDirtyKmers = true;
                    }
                }

                boolean isAlignedSomewhere = false;
                for (Interval iv : finalPos) {
                    if (!iv.getSequence().equals("*")) {
                        isAlignedSomewhere = true;
                    }
                }

                log.info("    alignment:");
                log.info("    - {}", finalPos);

                // See how many novel kmers we've used up
                int novelKmersContained = 0;
                int novelKmersUsed = 0;
                boolean hasRejectedKmers = false;
                boolean noConnections = false;

                for (AnnotatedVertex av : ag.vertexSet()) {
                    CortexKmer ck = new CortexKmer(av.getKmer());

                    if (novelKmers.containsKey(ck)) {
                        novelKmersContained++;

                        if (novelKmers.get(ck)) {
                            totalNovelKmersUsed++;
                            novelKmersUsed++;
                            novelKmers.put(ck, false);
                        }

                        CortexRecord nr = CLEAN.findRecord(ck);
                        if (nr == null) {
                            nr = DIRTY.findRecord(ck);
                        }

                        if (nr != null) {
                            if (nr.getInEdgesAsStrings(0).size() == 0 || nr.getOutEdgesAsStrings(0).size() == 0) {
                                noConnections = true;
                            }
                        }
                    }

                    if (rejects.contains(ck)) {
                        hasRejectedKmers = true;
                    }
                }

                if (ag.vertexSet().size() == 0) {
                    for (int i = 0; i <= stretch.length() - CLEAN.getKmerSize(); i++) {
                        String sk = stretch.substring(i, i + CLEAN.getKmerSize());
                        CortexKmer ck = new CortexKmer(sk);

                        if (novelKmers.containsKey(ck)) {
                            novelKmersContained++;

                            if (novelKmers.get(ck)) {
                                totalNovelKmersUsed++;
                                novelKmersUsed++;
                                novelKmers.put(ck, false);
                            }
                        }
                    }
                }

                gvc.attribute(0, "novelKmersContained", novelKmersContained);
                gvc.attribute(0, "novelKmersUsed", novelKmersUsed);

                log.info("    novelty:");
                log.info("    - novel kmers present:  {}/{}", novelKmersContained, novelKmers.size());
                log.info("    - novel kmers utilized: {}/{}", novelKmersUsed, novelKmers.size());
                log.info("    - cumulative usage:     {}/{}", totalNovelKmersUsed, novelKmers.size());

                //gvc.attribute(0, "filter", (novelKmersContained <= 1 && gvc.getAttributeAsString(0, "event").equals("unknown")) ? "FAIL" : "PASS");
                gvc.attribute(0, "filter", "PASS");

                String eventType = gvc.getAttributeAsString(0, "event");

                if (eventType.equals("unknown") || eventType.equals("GC") || eventType.equals("NAHR")) {
                    if (hasRejectedKmers) {
                        gvc.attribute(0, "filter", "HAS_REJECTED_KMERS");
                    } else if (!isAlignedSomewhere) {
                        gvc.attribute(0, "filter", "NO_ALIGNMENTS");
                    } else if (hasDirtyKmers) {
                        gvc.attribute(0, "filter", "HAS_DIRTY_KMERS");
                    } else if (noConnections) {
                        gvc.attribute(0, "filter", "DISCONNECTED");
                    }
                }

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

                Set<GraphicalVariantContext> newGvcs = new LinkedHashSet<GraphicalVariantContext>();

                // Separate MNPs into separate events if necessary
                if (gvc.getAttribute(0, "event").equals("MNP") && gvc.getAttributeAsString(0, "childAllele").length() == gvc.getAttributeAsString(0, "parentalAllele").length()) {
                    String parentalAllele = gvc.getAttributeAsString(0, "parentalAllele");
                    String childAllele = gvc.getAttributeAsString(0, "childAllele");

                    boolean[] mma = new boolean[parentalAllele.length()];
                    for (int i = 0; i < parentalAllele.length(); i++) {
                        mma[i] = (parentalAllele.charAt(i) != childAllele.charAt(i));
                    }

                    int startPos = 0;
                    StringBuilder sap = new StringBuilder();
                    StringBuilder sac = new StringBuilder();
                    for (int i = 0; i < parentalAllele.length(); i++) {
                        if (mma[i]) {
                            sap.append(parentalAllele.charAt(i));
                            sac.append(childAllele.charAt(i));
                        } else {
                            if (sap.length() > 0) {
                                GraphicalVariantContext newGvc = new GraphicalVariantContext(gvc);
                                newGvc.attribute(0, "event", sap.length() > 1 ? "MNP" : "SNP");
                                newGvc.attribute(0, "parentalAllele", sap.toString());
                                newGvc.attribute(0, "childAllele", sac.toString());
                                newGvcs.add(newGvc);

                                log.info("{} {} {}", startPos, sap.toString(), sac.toString());

                                sap = new StringBuilder();
                                sac = new StringBuilder();
                            }
                            startPos = i + 1;
                        }
                    }

                    if (sap.length() > 0) {
                        GraphicalVariantContext newGvc = new GraphicalVariantContext(gvc);
                        newGvc.attribute(0, "event", sap.length() > 1 ? "MNP" : "SNP");
                        newGvc.attribute(0, "parentalAllele", sap.toString());
                        newGvc.attribute(0, "childAllele", sac.toString());
                        newGvcs.add(newGvc);

                        log.info("{} {} {}", startPos, sap.toString(), sac.toString());
                    }
                } else {
                    newGvcs.add(gvc);
                }

                gvcs.addAll(newGvcs);

                int id = 0;
                for (GraphicalVariantContext newGvc : newGvcs) {
                    String key = novelKmer.getKmerAsString() + id;
                    evalTables.getTable("variantCalls").set(key, "novelKmer", novelKmer);
                    evalTables.getTable("variantCalls").set(key, "stretchNum", stretchNum);
                    evalTables.getTable("variantCalls").set(key, "filter", newGvc.getAttributeAsString(0, "filter"));
                    evalTables.getTable("variantCalls").set(key, "event", newGvc.getAttributeAsString(0, "event"));
                    evalTables.getTable("variantCalls").set(key, "isPolymorphic", newGvc.getAttributeAsBoolean(0, "isPolymorphic"));
                    evalTables.getTable("variantCalls").set(key, "locus", compactAlignments(finalPos));
                    evalTables.getTable("variantCalls").set(key, "parentalAlleleLength", newGvc.getAttributeAsString(0, "parentalAllele").length());
                    evalTables.getTable("variantCalls").set(key, "childAlleleLength", newGvc.getAttributeAsString(0, "childAllele").length());
                    evalTables.getTable("variantCalls").set(key, "parentalAllele", newGvc.getAttributeAsString(0, "parentalAllele"));
                    evalTables.getTable("variantCalls").set(key, "childAllele", newGvc.getAttributeAsString(0, "childAllele"));
                    evalTables.getTable("variantCalls").set(key, "score", newGvc.getAttributeAsInt(0, "score"));
                    evalTables.getTable("variantCalls").set(key, "haplotypeBackground", newGvc.getAttributeAsInt(0, "haplotypeBackground"));
                    evalTables.getTable("variantCalls").set(key, "start", newGvc.getAttributeAsInt(0, "start"));
                    evalTables.getTable("variantCalls").set(key, "stop", newGvc.getAttributeAsInt(0, "stop"));
                    evalTables.getTable("variantCalls").set(key, "novelKmersContained", novelKmersContained);
                    evalTables.getTable("variantCalls").set(key, "novelKmersUsed", novelKmersUsed);
                    evalTables.getTable("variantCalls").set(key, "novelKmersTotal", novelKmers.size());
                    evalTables.getTable("variantCalls").set(key, "traversalStatus", newGvc.getAttributeAsString(0, "traversalStatus"));
                    evalTables.getTable("variantCalls").set(key, "stretch", stretch);
                    evalTables.getTable("variantCalls").set(key, "parentalStretch", newGvc.getAttributeAsString(0, "parentalStretch"));
                    evalTables.getTable("variantCalls").set(key, "childStretch", newGvc.getAttributeAsString(0, "childStretch"));

                    if (finalPos.size() == 1) {
                        Interval pos = finalPos.get(0);
                        int radius = gvc.getAttributeAsInt(0, "haplotypeBackground") <= 1 ? 8 : 6;

                        if (gvc.getAttribute(0, "event").equals("GC")) {
                            cout.printf("stretch%d %s %d %d radius1=0.%dr # %s_%s\n", stretchNum, pos.getSequence(), pos.getStart(), pos.getEnd(), 8, gvc.getAttribute(0, "event"), gvc.getAttribute(0, "traversalStatus"));
                            cout.printf("stretch%d %s %d %d radius2=0.%dr # %s_%s\n", stretchNum, pos.getSequence(), pos.getStart(), pos.getEnd(), 6, gvc.getAttribute(0, "event"), gvc.getAttribute(0, "traversalStatus"));
                        } else {
                            cout.printf("stretch%d %s %d %d radius1=0.%dr # %s_%s\n", stretchNum, pos.getSequence(), pos.getStart(), pos.getEnd(), radius, gvc.getAttribute(0, "event"), gvc.getAttribute(0, "traversalStatus"));
                            cout.printf("stretch%d %s %d %d radius2=0.%dr # %s_%s\n", stretchNum, pos.getSequence(), pos.getStart(), pos.getEnd(), radius, gvc.getAttribute(0, "event"), gvc.getAttribute(0, "traversalStatus"));
                        }
                    } else if (finalPos.size() > 1) {
                        for (int i = 0; i < finalPos.size() - 1; i++) {
                            Interval sr1 = finalPos.get(i);
                            Interval sr2 = finalPos.get(i + 1);

                            int r1 = sr1.getName().equals(".") || Integer.valueOf(sr1.getName()) == 1 ? 8 : 6;
                            int r2 = sr2.getName().equals(".") || Integer.valueOf(sr2.getName()) == 1 ? 8 : 6;

                            cout.printf("stretch%d.%d %s %d %d radius1=0.%dr # %s_%s\n", stretchNum, i, sr1.getSequence(), sr1.getStart(), sr1.getEnd(), r1, gvc.getAttribute(0, "event"), gvc.getAttribute(0, "traversalStatus"));
                            cout.printf("stretch%d.%d %s %d %d radius2=0.%dr # %s_%s\n", stretchNum, i, sr2.getSequence(), sr2.getStart(), sr2.getEnd(), r2, gvc.getAttribute(0, "event"), gvc.getAttribute(0, "traversalStatus"));
                        }
                    } else {
                        cout.printf("stretch%d %s %d %d radius1=0.%dr # %s_%s\n", stretchNum, "NA", 100*stretchNum, 100*stretchNum, 7, gvc.getAttribute(0, "event"), gvc.getAttribute(0, "traversalStatus"));
                        cout.printf("stretch%d %s %d %d radius2=0.%dr # %s_%s\n", stretchNum, "NA", 100*stretchNum, 100*stretchNum, 7, gvc.getAttribute(0, "event"), gvc.getAttribute(0, "traversalStatus"));
                    }
                    id++;
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
