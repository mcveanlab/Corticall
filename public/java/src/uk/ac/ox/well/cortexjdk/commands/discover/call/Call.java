package uk.ac.ox.well.cortexjdk.commands.discover.call;

import com.google.common.base.Joiner;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;
import org.jgrapht.Graphs;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.mosaic.MosaicAligner;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.alignment.sw.SmithWaterman;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.containers.ContainerUtils;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContigStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.DestinationStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.*;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

import static htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder.OutputType.VCF;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.REVERSE;

public class Call extends Module {
    @Argument(fullName = "graph", shortName = "g", doc = "Graph")
    public CortexGraph GRAPH;

    @Argument(fullName = "links", shortName = "l", doc = "Links", required=false)
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName = "rois", shortName = "r", doc = "Rois")
    public CortexGraph ROIS;

    @Argument(fullName="partitions", shortName="p", doc="Partitions")
    public FastaSequenceFile PARTITIONS;

    @Argument(fullName="mother", shortName="m", doc="Mother's sample name")
    public LinkedHashSet<String> MOTHER;

    @Argument(fullName="father", shortName="f", doc="Father's sample name")
    public LinkedHashSet<String> FATHER;

    @Argument(fullName="background", shortName="b", doc="Background", required=false)
    public HashMap<String, IndexedReference> BACKGROUNDS;

    @Argument(fullName="partitionName", shortName="pn", doc="Partitions to process", required=false)
    public HashSet<String> PARTITION_NAMES;

    //@Argument(fullName="reference", shortName="R", doc="Reference", required=false)
    //public HashMap<String, IndexedReference> REFERENCE;

    @Output
    public File out;

    @Output(fullName="accountingOut", shortName="ao", doc="Accounting out")
    public PrintStream aout;

    @Override
    public void execute() {
        Set<CanonicalKmer> rois = loadRois(ROIS);
        List<ReferenceSequence> rseqs = loadPartitions();

        SAMSequenceDictionary sd = buildMergedSequenceDictionary(rseqs);
        Set<VariantContext> svcs = buildVariantSorter(sd);
        VariantContextWriter vcw = buildVariantWriter(sd);

        MosaicAligner ma = new MosaicAligner(0.35, 0.99, 1e-5, 0.001);

        for (int rseqIndex = 0; rseqIndex < rseqs.size(); rseqIndex++) {
            ReferenceSequence rseq = rseqs.get(rseqIndex);
            String seq = rseq.getBaseString();

            List<CortexVertex> w = loadChildWalk(rseq, GRAPH);
            List<Triple<Integer, Integer, List<CortexVertex>>> sections = sectionContig(rois, w, 100, 500);

            Set<VariantContextBuilder> vcs = buildVariantContextBuilderSorter(sd);

            if (sections == null) {
                log.info("Partition {}/{} (sections={}, fullname={}) [skipped]", rseqIndex, rseqs.size() - 1, 0, rseq.getName());
            } else {
                log.info("Partition {}/{} (sections={}, fullname={})", rseqIndex, rseqs.size() - 1, sections.size(), rseq.getName());

                for (int sectionIndex = 0; sectionIndex < sections.size(); sectionIndex++) {
                    log.debug("  section {}/{}", sectionIndex + 1, sections.size());

                    Triple<Integer, Integer, List<CortexVertex>> section = sections.get(sectionIndex);

                    List<CortexVertex> ws = section.getRight();

                    Map<String, Map<String, String>> allTargets = new HashMap<>();
                    allTargets.put("all", new HashMap<>());
                    for (Set<String> parentName : Arrays.asList(MOTHER, FATHER)) {
                        Map<String, String> parentalTargets = fasterAssembleCandidateHaplotypes(ws, parentName);

                        allTargets.get("all").putAll(parentalTargets);
                        allTargets.put(parentName.iterator().next(), parentalTargets);
                    }

                    if (allTargets.get("all").size() > 0) {
                        Triple<Integer, Integer, String> trimmedQuery = trimQuery(ws, allTargets.get("all"), rois);

                        String bestTagName = "all";
                        double llk = Double.MIN_VALUE;
                        for (String tagName : allTargets.keySet()) {
                            if (allTargets.get(tagName).size() > 0) {
                                ma.align(trimmedQuery.getRight(), allTargets.get(tagName));

                                if (ma.getMaximumLogLikelihood() > llk) {
                                    bestTagName = tagName;
                                    llk = ma.getMaximumLogLikelihood();
                                }
                            }
                        }

                        Map<String, String> targets = allTargets.get(bestTagName);

                        List<Triple<String, String, Pair<Integer, Integer>>> lps = ma.align(trimmedQuery.getRight(), targets);
                        log.debug("\n{}\n{}", makeNoveltyTrack(rois, lps, true), ma);

                        List<Pair<Integer, Integer>> nrs = getNoveltyRegions(rois, lps, true);

                        List<VariantContextBuilder> calls = new ArrayList<>();
                        calls.addAll(callSmallBubbles(lps, nrs, rseq.getName().split(" ")[0], section.getLeft() + trimmedQuery.getLeft(), section.getMiddle() + trimmedQuery.getLeft()));
                        calls.addAll(callLargeBubbles(lps, nrs, targets, rseq.getName().split(" ")[0], section.getLeft() + trimmedQuery.getLeft(), section.getMiddle() + trimmedQuery.getLeft()));
                        calls.addAll(callBreakpoints(lps, nrs, rseq.getName().split(" ")[0], section.getLeft() + trimmedQuery.getLeft(), section.getMiddle() + trimmedQuery.getLeft()));

                        List<VariantContextBuilder> merged = mergeBubbles(lps, calls);

                        Set<CanonicalKmer> sectionRois = new TreeSet<>();
                        for (int i = 0; i <= trimmedQuery.getRight().length() - GRAPH.getKmerSize(); i++) {
                            CanonicalKmer ck = new CanonicalKmer(trimmedQuery.getRight().substring(i, i + GRAPH.getKmerSize()));
                            if (rois.contains(ck)) {
                                sectionRois.add(ck);
                            }
                        }

                        for (VariantContextBuilder vcb : merged) {
                            vcb.attribute("targets", targets);
                            vcb.attribute("lps", lps);
                            vcb.attribute("sectionIndex", sectionIndex);
                            vcb.attribute("novels", Joiner.on(",").join(sectionRois));
                        }

                        vcs.addAll(merged);
                    }
                }
            }

            vcs = filterBreakpoints(vcs);

            vcs = mergeBreakpoints(seq, vcs, rois);

            vcs = mergeDoubleBreakpoints(seq, vcs);
            //vcs = mergeSingleBreakpoints(seq, vcs);

            vcs = assignCoordinates(vcs);

            for (VariantContextBuilder vcb : vcs) {
                vcb.rmAttributes(Arrays.asList(
                        "targets", "lps",
                        "nextBase", "nextChrom", "nextStart", "nextStop", "nextStrand",
                        "prevBase", "prevChrom", "prevStart", "prevStop", "prevStrand",
                        "targetName", "targetStart", "targetStop",
                        "start", "stop",
                        "sectionStart", "sectionStop",
                        "variantStart", "variantStop"
                ));

                VariantContext vc = vcb.make();

                if (!vc.isFiltered()) {
                    String back = vc.getAttributeAsString("BACKGROUND", "unknown");
                    int start = vc.getStart();
                    int end = vc.isSymbolic() ? start : vc.getEnd();
                    String refAllele = BACKGROUNDS.containsKey(back) && BACKGROUNDS.get(back).getReferenceSequence().getSequenceDictionary().getSequence(vc.getContig()) != null ? BACKGROUNDS.get(back).getReferenceSequence().getSubsequenceAt(vc.getContig(), start, end).getBaseString() : "?";
                    log.debug("{} {} {} {}", vc.getReference(), refAllele, vc.getFilters(), new VariantContextBuilder(vc).rmAttribute("novels").make());

                    svcs.add(vc);
                }
            }

            log.debug("");
        }

        writeVariants(rois, svcs, vcw);
    }

    @NotNull
    private Set<VariantContextBuilder> filterBreakpoints(Set<VariantContextBuilder> vcs) {
        IntervalTreeMap<List<VariantContextBuilder>> itm = new IntervalTreeMap<>();

        for (VariantContextBuilder vcbn : vcs) {
            VariantContext vc = vcbn.make();

            Interval vit = new Interval(vc.getContig(), vc.getStart(), vc.getEnd());
            if (!itm.containsKey(vit)) {
                itm.put(vit, new ArrayList<>());
            }
            itm.get(vit).add(vcbn);
        }

        Set<String> toFilter = new HashSet<>();
        for (Interval it : itm.keySet()) {
            boolean hasConcreteVariants = false, hasBreakpointVariants = false;
            if (itm.get(it).size() > 1) {
                for (VariantContextBuilder vc : itm.get(it)) {
                    if (!vc.make().isSymbolic()) {
                        hasConcreteVariants = true;
                    } else if (vc.make().hasAttribute("MATEID")) {
                        hasBreakpointVariants = true;
                    }
                }
            }

            if (hasConcreteVariants && hasBreakpointVariants) {
                for (VariantContextBuilder vc : itm.get(it)) {
                    if (vc.make().isSymbolic() && vc.make().hasAttribute("MATEID")) {
                        toFilter.add(vc.make().getID());
                        toFilter.add(vc.make().getAttributeAsString("MATEID", "unknown"));
                    }
                }
            }
        }

        for (Interval it : itm.keySet()) {
            for (VariantContextBuilder vc : itm.get(it)) {
                if (toFilter.contains(vc.make().getID())) {
                    vc.filter("OVERLAPPING_BREAKPOINT");
                }
            }
        }

        Set<VariantContextBuilder> newVcbs = new HashSet<>();
        for (List<VariantContextBuilder> vcbs : itm.values()) {
            newVcbs.addAll(vcbs);
        }

        return newVcbs;
    }

    private Set<VariantContextBuilder> assignCoordinates(Set<VariantContextBuilder> vcs) {
        Set<VariantContextBuilder> newVcbs = new HashSet<>();
        for (VariantContextBuilder vcb : vcs) {
            newVcbs.add(assignCoordinates(vcb));
        }

        return newVcbs;
    }

    private VariantContextBuilder assignCoordinates(VariantContextBuilder vcb) {
        VariantContextBuilder vcbn = new VariantContextBuilder(vcb);

        List<Triple<String, String, Pair<Integer, Integer>>> lps = (ArrayList<Triple<String, String, Pair<Integer, Integer>>>) vcbn.make().getAttribute("lps");

        int start = vcbn.make().getAttributeAsInt("start", 0) + (vcbn.make().isSNP() ? 1 : 0);
        StringBuilder prevFlankBuilder = new StringBuilder();
        for (int q = start; q >= 0 && getParentalRow(lps, start) == getParentalRow(lps, q); q--) {
            if (getParentalColumn(lps, q) != '-') {
                prevFlankBuilder.insert(0, getParentalColumn(lps, q));
            }
        }

        String prevBack = lps.get(getParentalRow(lps, start)).getLeft().split(":")[0];
        String prevFlank = prevFlankBuilder.toString();
        List<SAMRecord> prevSrs = sortAlignments(prevBack, prevFlank);
        SAMRecord prevSr = prevSrs.size() == 0 ? null : prevSrs.get(0);
        if (prevSr != null) {
            vcbn.attribute("prevChrom", prevSr.getContig());
            vcbn.attribute("prevStart", prevSr.getReferencePositionAtReadPosition(1) + 1);
            vcbn.attribute("prevStop", prevSr.getReferencePositionAtReadPosition(prevSr.getReadLength()) + 1);
            vcbn.attribute("prevStrand", prevSr.getReadNegativeStrandFlag() ? "-" : "+");
        }

        int stop = vcbn.make().getAttributeAsInt("stop", 0) - (vcbn.make().isSNP() ? 1 : 0);
        while (getParentalColumn(lps, stop) == '-' && stop < lps.get(0).getMiddle().length()) {
            stop++;
        }
        StringBuilder nextFlankBuilder = new StringBuilder();
        for (int q = stop; q < lps.get(0).getMiddle().length() && getParentalRow(lps, stop) == getParentalRow(lps, q); q++) {
            if (getParentalColumn(lps, q) != '-') {
                nextFlankBuilder.append(getParentalColumn(lps, q));
            }
        }

        String nextBack = lps.get(getParentalRow(lps, stop)).getLeft().split(":")[0];
        String nextFlankFw = nextFlankBuilder.toString();
        List<SAMRecord> nextSrs = sortAlignments(nextBack, nextFlankFw);
        SAMRecord nextSr = nextSrs.size() == 0 ? null : nextSrs.get(0);

        if (prevSr != null && nextSrs.size() > 0) {
            for (SAMRecord nsr : nextSrs) {
                if (prevSr.getContig().equals(nsr.getContig())) {
                    nextSr = nsr;
                    break;
                }
            }
        }

        if (nextSr != null) {
            vcbn.attribute("nextChrom", nextSr.getContig());
            vcbn.attribute("nextStart", nextSr.getReferencePositionAtReadPosition(1) + 1);
            vcbn.attribute("nextStop", nextSr.getReferencePositionAtReadPosition(nextSr.getReadLength()) + 1);
            vcbn.attribute("nextStrand", nextSr.getReadNegativeStrandFlag() ? "-" : "+");
        }

        SAMRecord sr = null;
        List<SAMRecord> srs = null;
        int alignStart = 0;
        if (prevSr != null && nextSr != null) {
            //sr = prevSr.getStart() < nextSr.getStart() ? prevSr : nextSr;
            //srs = prevSr.getStart() < nextSr.getStart() ? prevSrs : nextSrs;
            //alignStart = sr.getEnd() + 1;

            if (prevSr.getStart() < nextSr.getStart()) {
                nextSr = null;
            } else {
                prevSr = null;
            }
        }

        if (prevSr != null) {
            sr = prevSr;
            srs = prevSrs;
            alignStart = sr.getReadNegativeStrandFlag() ? sr.getStart() + 1 : sr.getEnd() + 1;
        } else if (nextSr != null) {
            sr = nextSr;
            srs = nextSrs;
            alignStart = sr.getReadNegativeStrandFlag() ? sr.getEnd() + 1 : sr.getStart();
        }

        if (sr != null) {
            boolean flip = sr.getReadNegativeStrandFlag();
            List<Allele> alleles = vcbn.getAlleles();

            vcbn.chr(sr.getContig());
            vcbn.start(alignStart);
            vcbn.stop(alignStart + vcb.make().getEnd() - vcb.make().getStart());
            vcbn.attribute("flankMappingQuality", sr.getMappingQuality());

            if (flip) {
                List<Allele> allelesRc = new ArrayList<>();

                for (Allele a : alleles) {
                    String[] pieces = new String(a.getDisplayBases()).split("((?<=[\\[\\]])|(?=[\\[\\]]))");
                    for (int i = 0; i < pieces.length; i++) {
                        String piece = pieces[i];

                        if (piece.matches("^(\\.?)[ACTGacgt]+(\\.?)$")) {
                            piece = SequenceUtils.reverseComplement(piece);
                            pieces[i] = piece;
                        }
                    }

                    String newAllele = Joiner.on("").join(pieces);

                    if (!vcbn.make().isSNP() && !vcbn.make().isSymbolic()) {
                        String newRefBase = SequenceUtils.reverseComplement(sr.getReadString().substring(0, 1));
                        newAllele = newRefBase + newAllele.substring(0, newAllele.length() - 1);
                    }

                    allelesRc.add(Allele.create(newAllele, a.isReference()));
                }

                alleles = allelesRc;
            }

            List<Allele> allelesRevised = new ArrayList<>();
            for (int i = 0; i < alleles.size(); i++) {
                Allele a = alleles.get(i);
                String[] pieces = new String(a.getDisplayBases()).split("((?<=[\\[\\]])|(?=[\\[\\]]))");

                if (pieces.length == 4) {
                    String[] newpieces = new String[4];

                    if (pieces[3].matches("^(\\.?)[ACTGacgt]+(\\.?)$")) {
                        newpieces[0] = pieces[3];
                        newpieces[1] = pieces[0].equals("[") ? "]" : "[";
                        newpieces[2] = pieces[1];
                        newpieces[3] = pieces[2].equals("[") ? "]" : "[";

                        String[] subpieces = pieces[1].split(":");
                        String back = subpieces[0];
                        String contigName = subpieces[0] + ":" + subpieces[1] + ":" + subpieces[2];

                        for (int m = 1; m < lps.size(); m++) {
                            if (lps.get(m).getLeft().equals(contigName)) {
                                if (BACKGROUNDS.containsKey(back)) {
                                    List<SAMRecord> mrs = sortAlignments(back, lps.get(m).getMiddle().replaceAll(" ", ""));

                                    if (mrs.size() > 0) {
                                        SAMRecord mr = mrs.get(0);
                                        int newpos = mr.getReferencePositionAtReadPosition(1);

                                        newpieces[2] = mr.getContig() + ":" + newpos;
                                    }
                                }

                                break;
                            }
                        }
                    } else {
                        newpieces[0] = pieces[1].equals("[") ? "]" : "[";
                        newpieces[1] = pieces[2];
                        newpieces[2] = pieces[3].equals("[") ? "]" : "[";
                        newpieces[3] = pieces[0];

                        String[] subpieces = pieces[2].split(":");
                        String back = subpieces[0];
                        String contigName = subpieces[0] + ":" + subpieces[1] + ":" + subpieces[2];

                        for (int m = 1; m < lps.size(); m++) {
                            if (lps.get(m).getLeft().equals(contigName)) {
                                if (BACKGROUNDS.containsKey(back)) {
                                    List<SAMRecord> mrs = sortAlignments(back, lps.get(m).getMiddle().replaceAll(" ", ""));

                                    if (mrs.size() > 0) {
                                        SAMRecord mr = mrs.get(0);
                                        int newpos = mr.getReferencePositionAtReadPosition(1);

                                        newpieces[1] = mr.getContig() + ":" + newpos;
                                    }
                                }

                                break;
                            }
                        }
                    }

                    pieces = newpieces;

                    allelesRevised.add(Allele.create(Joiner.on("").join(pieces)));
                } else {
                    allelesRevised.add(a);
                }
            }

            vcbn.alleles(allelesRevised);
            vcbn.attribute("flipped", flip);

            List<String> altLoci = new ArrayList<>();
            for (SAMRecord sra : srs) {
                String locus = String.format("%s:%d", sra.getContig(), sra.getStart(), sra.getReadNegativeStrandFlag());

                altLoci.add(locus);
            }

            vcbn.attribute("alt_loci", Joiner.on(",").join(altLoci));
        }

        return vcbn;
    }

    private List<VariantContextBuilder> mergeBubbles(List<Triple<String, String, Pair<Integer, Integer>>> lps, List<VariantContextBuilder> calls) {
        if (calls.size() <= 1) {
            return calls;
        }

        List<VariantContextBuilder> merged = new ArrayList<>();

        for (int i = 0; i < calls.size(); i++) {
            if (i + 1 <= calls.size() - 1) {
                int start0 = calls.get(i).make().getAttributeAsInt("start", 0);
                int stop0 = calls.get(i).make().getAttributeAsInt("stop", 0);
                int start1 = calls.get(i + 1).make().getAttributeAsInt("start", 500);
                int stop1 = calls.get(i + 1).make().getAttributeAsInt("stop", 500);

                if (start1 - stop0 < 10 && !calls.get(i).make().isSymbolicOrSV() && !calls.get(i+1).make().isSymbolicOrSV()) {
                    StringBuilder cBuilder = new StringBuilder();
                    StringBuilder pBuilder = new StringBuilder();

                    for (int j = start0; j < stop1; j++) {
                        char c = getChildColumn(lps, j);
                        char p = getParentalColumn(lps, j);

                        if (c != '-') {
                            cBuilder.append(c);
                        }
                        if (p != '-') {
                            pBuilder.append(p);
                        }
                    }

                    List<Allele> alleles = Arrays.asList(Allele.create(pBuilder.toString(), true), Allele.create(cBuilder.toString()));

                    char prevBase = getChildColumn(lps, start0);
                    char nextBase = getChildColumn(lps, stop1);

                    int sectionStart = calls.get(i).make().getAttributeAsInt("sectionStart", 0);

                    VariantContextBuilder vcb = new VariantContextBuilder(calls.get(i))
                            .alleles(alleles)
                            .start(sectionStart + start0)
                            .computeEndFromAlleles(alleles, sectionStart + start0)
                            .attribute("start", start0)
                            .attribute("stop", stop1)
                            .attribute("variantStart", sectionStart + start0)
                            .attribute("variantStop", sectionStart + stop1)
                            .attribute("prevBase", prevBase)
                            .attribute("nextBase", nextBase);

                    if (cBuilder.substring(1, cBuilder.length()).equals(SequenceUtils.reverseComplement(pBuilder.substring(1, pBuilder.length())))) {
                        vcb.attribute("SVTYPE", "INV");
                    }

                    merged.add(vcb);

                    i = i + 1;
                } else {
                    merged.add(calls.get(i));
                }
            } else {
                merged.add(calls.get(i));
            }
        }

        return merged;
    }

    private Set<VariantContextBuilder> mergeSingleBreakpoints(String seq, Set<VariantContextBuilder> callset) {
        List<VariantContextBuilder> calls = new ArrayList<>(callset);

        if (calls.size() <= 1) {
            return callset;
        }

        List<VariantContextBuilder> bnds = new ArrayList<>();

        for (int i = 0; i < calls.size(); i++) {
            VariantContext vc = calls.get(i).make();
            if (vc.isSymbolicOrSV() && vc.getAttributeAsString("SVTYPE", "unknown").equals("BND")) {
                bnds.add(calls.get(i));
            }
        }

        Map<String, VariantContextBuilder> replacements = new HashMap<>();
        Set<Interval> removals = new HashSet<>();

        for (int i = 0; i < bnds.size() - 1; i++) {
            VariantContextBuilder vcb0 = bnds.get(i);
            VariantContext vc0 = vcb0.make();

            String back = vc0.getAttributeAsString("BACKGROUND", "");

            for (int j = i + 1; j < bnds.size(); j++) {
                VariantContextBuilder vcb1 = bnds.get(j);
                VariantContext vc1 = vcb1.make();

                List<Triple<String, String, Pair<Integer, Integer>>> lps0 = (ArrayList<Triple<String, String, Pair<Integer, Integer>>>) vc0.getAttribute("lps");
                int q0 = 0;
                StringBuilder kmer0 = new StringBuilder();
                for (int q = 1; q < lps0.size(); q++) {
                    if (lps0.get(q).getLeft().equals(vc0.getAttributeAsString("targetName", ""))) {
                        int pos = vc0.getAttributeAsInt("start", 0);
                        while (kmer0.length() < GRAPH.getKmerSize() && pos >= 0) {
                            char c = lps0.get(q).getMiddle().charAt(pos);

                            if (c != '-') {
                                kmer0.insert(0, c);
                            }

                            pos--;
                        }

                        q0 = q;

                        break;
                    }
                }

                List<Triple<String, String, Pair<Integer, Integer>>> lps1 = (ArrayList<Triple<String, String, Pair<Integer, Integer>>>) vc1.getAttribute("lps");
                int q1 = 0;
                StringBuilder kmer1 = new StringBuilder();
                for (int q = 1; q < lps1.size(); q++) {
                    if (lps1.get(q).getLeft().equals(vc1.getAttributeAsString("targetName", ""))) {
                        int pos = vc1.getAttributeAsInt("start", 0);
                        while (kmer1.length() < GRAPH.getKmerSize() && pos <= lps1.get(q).getMiddle().length()) {
                            char c = lps1.get(q).getMiddle().charAt(pos);

                            if (c != '-') {
                                kmer1.append(c);
                            }

                            pos++;
                        }

                        q1 = q;

                        break;
                    }
                }

                for (Set<String> parentName : Arrays.asList(MOTHER, FATHER)) {
                    TraversalEngine ef = new TraversalEngineFactory()
                            .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                            .traversalDirection(FORWARD)
                            .combinationOperator(OR)
                            .stoppingRule(DestinationStopper.class)
                            .graph(GRAPH)
                            .links(LINKS)
                            .make();

                    TraversalEngine er = new TraversalEngineFactory()
                            .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                            .traversalDirection(REVERSE)
                            .combinationOperator(OR)
                            .stoppingRule(DestinationStopper.class)
                            .graph(GRAPH)
                            .links(LINKS)
                            .make();

                    DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = ef.dfs(kmer0.toString(), kmer1.toString());
                    if (g == null || g.vertexSet().size() == 0) {
                        g = er.dfs(kmer1.toString(), kmer0.toString());
                    }
                }

                if (BACKGROUNDS.containsKey(back) && vc0.getContig().equals(vc1.getContig()) && vc0.getAttributeAsString("MATEID", "").equals(vc1.getID())) {
                    int pStart = vc0.getStart() < vc1.getStart() ? vc0.getStart() : vc1.getStart();
                    int pEnd = vc0.getEnd() > vc1.getEnd() ? vc0.getEnd() : vc1.getEnd();
                    int vStart = vc0.getAttributeAsInt("variantStart", 0) < vc1.getAttributeAsInt("variantStart", 0) ? vc0.getAttributeAsInt("variantStart", 0) : vc1.getAttributeAsInt("variantStart", 0);
                    int vEnd = vc0.getAttributeAsInt("variantStop", 0) > vc1.getAttributeAsInt("variantStop", 0) ? vc0.getAttributeAsInt("variantStop", 0) : vc1.getAttributeAsInt("variantStop", 0);

                    String parentalContig = BACKGROUNDS.get(back).getReferenceSequence().getSubsequenceAt(vc0.getContig(), pStart, pEnd).getBaseString();
                    String childContig = seq.substring(vStart, vEnd);

                    if (vc0.getAttributeAsBoolean("flipped", false)) {
                        childContig = SequenceUtils.reverseComplement(childContig);
                    }

                    List<Allele> alleles = Arrays.asList(Allele.create(parentalContig, true), Allele.create(childContig));

                    VariantContextBuilder vcb = new VariantContextBuilder(vc0)
                            .alleles(alleles)
                            .start(pStart)
                            .computeEndFromAlleles(alleles, pStart)
                            .attribute("prevBase", vc0.getAttributeAsString("prevBase", "N"))
                            .attribute("nextBase", vc1.getAttributeAsString("nextBase", "N"))
                            .rmAttribute("SVTYPE")
                            .rmAttribute("MATEID");

                    replacements.put(vc0.getID(), vcb);
                    replacements.put(vc1.getID(), null);

                    removals.add(new Interval(vc0.getContig(), vc0.getStart(), vc0.getStart()));
                    removals.add(new Interval(vc1.getContig(), vc1.getStart(), vc1.getStart()));
                }
            }
        }

        Set<VariantContextBuilder> newcalls = new LinkedHashSet<>();
        for (VariantContextBuilder vcb : calls) {
            if (!vcb.make().isSymbolic() && removals.contains(new Interval(vcb.make().getContig(), vcb.make().getStart(), vcb.make().getStart()))) {
                // ignore
            } else {
                if (!replacements.containsKey(vcb.make().getID())) {
                    newcalls.add(vcb);
                } else {
                    if (replacements.get(vcb.make().getID()) != null) {
                        newcalls.add(replacements.get(vcb.make().getID()));
                    }
                }
            }
        }

        return newcalls;
    }

    private Set<VariantContextBuilder> mergeBreakpoints(String seq, Set<VariantContextBuilder> callset, Set<CanonicalKmer> rois) {
        List<VariantContextBuilder> calls = new ArrayList<>(callset);

        if (calls.size() <= 1) {
            return callset;
        }

        List<VariantContextBuilder> bnds = new ArrayList<>();

        for (int i = 0; i < calls.size(); i++) {
            VariantContext vc = calls.get(i).make();
            if (vc.isSymbolicOrSV() && vc.getAttributeAsString("SVTYPE", "unknown").equals("BND")) {
                bnds.add(calls.get(i));
            }
        }

        bnds.sort((vb0, vb1) -> {
            VariantContext v0 = vb0.make();
            VariantContext v1 = vb1.make();

            if (v0.getStart() == v1.getStart()) { return 0;  }
            if (v0.getStart() <  v1.getStart()) { return -1; }
            return 1;
        });

        Set<VariantContextBuilder> removals = new HashSet<>();
        Set<VariantContextBuilder> additions = new HashSet<>();

        for (int i = 0; i < bnds.size() - 1; i += 2) {
            for (int j = i + 1; j < bnds.size(); j += 2) {
                VariantContextBuilder vb0 = bnds.get(i);
                VariantContext v0 = vb0.make();
                List<Triple<String, String, Pair<Integer, Integer>>> lps0 = (ArrayList<Triple<String, String, Pair<Integer, Integer>>>) v0.getAttribute("lps");

                VariantContextBuilder vb1 = bnds.get(j);
                VariantContext v1 = vb1.make();
                List<Triple<String, String, Pair<Integer, Integer>>> lps1 = (ArrayList<Triple<String, String, Pair<Integer, Integer>>>) v1.getAttribute("lps");

                String back0 = v0.getAttributeAsString("BACKGROUND", "");
                String back1 = v1.getAttributeAsString("BACKGROUND", "");

                String flank0 = lps0.get(0).getMiddle().substring(0, v0.getAttributeAsInt("start", 0) + 1).replaceAll("[- ]", "");
                String flank1 = lps1.get(0).getMiddle().substring(v1.getAttributeAsInt("start", 0), lps1.get(0).getMiddle().length()).replaceAll("[- ]", "");

                String kmer0 = flank0.length() - GRAPH.getKmerSize() - 1 >= 0 ? flank0.substring(flank0.length() - GRAPH.getKmerSize() - 1, flank0.length() - 1) : "";
                String kmer1 = GRAPH.getKmerSize() + 1 < flank1.length() ? flank1.substring(1, GRAPH.getKmerSize() + 1) : "";

                //String childContig = seq.substring(v0.getEnd() - GRAPH.getKmerSize() + 1, v1.getStart());
                int c0 = v0.getEnd() - GRAPH.getKmerSize() >= 0 ? v0.getEnd() - GRAPH.getKmerSize() : 0;
                int c1 = v1.getStart() + GRAPH.getKmerSize() + 1 < seq.length() ? v1.getStart() + GRAPH.getKmerSize() + 1 : seq.length();
                String childContig = seq.substring(c0, c1);
                String parentContig = null;

                log.info("{} {} {} {}", i, j, kmer0, kmer1);

                if (kmer0.length() == GRAPH.getKmerSize() && kmer1.length() == GRAPH.getKmerSize()) {
                    for (Set<String> parentName : Arrays.asList(MOTHER, FATHER)) {
                        if (parentName.contains(v0.getAttributeAsString("BACKGROUND", "unknown"))) {
                            TraversalEngine ef = new TraversalEngineFactory()
                                    .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                                    .traversalDirection(FORWARD)
                                    .combinationOperator(OR)
                                    .stoppingRule(DestinationStopper.class)
                                    .graph(GRAPH)
                                    .links(LINKS)
                                    .make();

                            TraversalEngine er = new TraversalEngineFactory()
                                    .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                                    .traversalDirection(REVERSE)
                                    .combinationOperator(OR)
                                    .stoppingRule(DestinationStopper.class)
                                    .graph(GRAPH)
                                    .links(LINKS)
                                    .make();

                            DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = ef.dfs(kmer0, kmer1);
                            if (g == null || g.vertexSet().size() == 0) {
                                g = er.dfs(kmer1, kmer0);
                            }

                            if (g != null && g.vertexSet().size() > 0) {
                                Set<String> parentalContigs = new HashSet<>();
                                ConnectivityInspector<CortexVertex, CortexEdge> ci = new ConnectivityInspector<>(g);
                                for (Set<CortexVertex> cvs : ci.connectedSets()) {
                                    List<CortexVertex> w = new ArrayList<>();

                                    for (CortexVertex cv : cvs) {
                                        List<CortexVertex> wa = TraversalUtils.toWalk(g, cv.getKmerAsString(), g.edgeSet().iterator().next().getColor());

                                        if (wa.size() == w.size()) {
                                            break;
                                        } else if (wa.size() > w.size()) {
                                            w = wa;
                                        }
                                    }

                                    if (w.size() > 0) {
                                        parentalContigs.add(TraversalUtils.toContig(w));
                                    }
                                }

                                for (String parentalContig : parentalContigs) {
                                    //log.info("{}", parentalContig);
                                    //log.info("{}", childContig);
                                    //log.info("");

                                    if (parentContig == null || parentalContig.length() > parentContig.length()) {
                                        parentContig = parentalContig;
                                    }
                                }
                            }
                        }
                    }
                }

                if (parentContig == null) {
                    List<SAMRecord> ss0 = sortAlignments(back0, flank0);
                    List<SAMRecord> ss1 = sortAlignments(back1, flank1);

                    SAMRecord s0 = ss0.size() > 0 ? ss0.get(0) : null;
                    SAMRecord s1 = ss1.size() > 0 ? ss1.get(0) : null;

                    if (s0 != null && s1 != null && back0.equals(back1) && s0.getContig().equals(s1.getContig()) && s0.getReadNegativeStrandFlag() == s1.getReadNegativeStrandFlag()) {
                        String contig = s0.getContig();
                        int left = Math.min(s1.getEnd() + 1, s0.getStart() - 1) + 1;
                        int right = Math.max(s1.getEnd() + 1, s0.getStart() - 1) + 1;

                        String parentalContig = BACKGROUNDS.get(back0).getReferenceSequence().getSubsequenceAt(contig, left, right).getBaseString();

                        parentContig = s0.getReadNegativeStrandFlag() ? SequenceUtils.reverseComplement(parentalContig) : parentalContig;
                    }
                }

                if (parentContig != null) {
                    float fedit = pctSharedKmers(childContig, parentContig);
                    float wedit = pctSharedKmers(SequenceUtils.reverse(childContig), parentContig);

                    //SmithWaterman sw = new SmithWaterman();
                    //String[] f = sw.getAlignment(childContig, parentContig);
                    //float fedit = 100.0f * ((float) f[0].length() - numEdits(f)) / (float) f[0].length();

                    //String[] w = sw.getAlignment(SequenceUtils.reverseComplement(childContig), parentContig);
                    //float wedit = 100.0f * ((float) w[0].length() - numEdits(w)) / (float) w[0].length();

                    //log.info("{}", fedit);
                    //log.info("{}", f[0]);
                    //log.info("{}", f[1]);

                    //log.info("{}", wedit);
                    //log.info("{}", w[0]);
                    //log.info("{}", w[1]);

                    if (fedit > wedit && fedit >= 50.0f && childContig.length() < 2000 && parentContig.length() < 2000) {
                        MosaicAligner ma = new MosaicAligner(0.35, 0.99, 1e-5, 0.001);
                        Map<String, String> newTargets = new HashMap<>();
                        newTargets.put(String.format("%s:%s_unknown:%s_contig0_merge", back0, back0, back0), parentContig);
                        List<Triple<String, String, Pair<Integer, Integer>>> newlps = ma.align(childContig, newTargets);

                        //List<Triple<String, String, Pair<Integer, Integer>>> newlps = new ArrayList<>();
                        //newlps.add(Triple.of("query", childContig, Pair.create(0, childContig.length())));
                        //newlps.add(Triple.of(String.format("%s:%s_unknown:%s_contig0_merge", back0, back0, back0), parentContig, Pair.create(0, childContig.length())));

                        List<Pair<Integer, Integer>> nrs = getNoveltyRegions(rois, newlps, true);

                        //log.info("\n{}\n{}", makeNoveltyTrack(rois, newlps, true), ma);

                        List<VariantContextBuilder> newCalls = new ArrayList<>();
                        newCalls.addAll(callSmallBubbles(newlps, nrs, v0.getAttributeAsString("PARTITION_NAME", ""), 0, childContig.length()));

                        newCalls = mergeBubbles(newlps, newCalls);

                        for (VariantContextBuilder vcb : newCalls) {
                            int ostart = vcb.make().getStart();
                            int ostop = vcb.make().getEnd();

                            additions.add(
                                    vcb
                                        .start(ostart + seq.indexOf(kmer0))
                                        .stop(ostop + seq.indexOf(kmer0))
                                        .attributes(v0.getAttributes())
                            );
                        }

                        //log.info("");
                        //calls.addAll(callLargeBubbles(newlps, nrs, targets, rseq.getName().split(" ")[0], section.getLeft() + trimmedQuery.getLeft(), section.getMiddle() + trimmedQuery.getLeft()));
                        //calls.addAll(callBreakpoints(newlps, nrs, rseq.getName().split(" ")[0], section.getLeft() + trimmedQuery.getLeft(), section.getMiddle() + trimmedQuery.getLeft()));
                    } else if (wedit > fedit && wedit >= 50.0f) {
                        List<Allele> alleles = Arrays.asList(Allele.create(parentContig, true), Allele.create(SequenceUtils.reverseComplement(childContig)));

                        //log.info("Inversion");

                        VariantContextBuilder vcb = new VariantContextBuilder(v0)
                                .alleles(alleles)
                                .computeEndFromAlleles(alleles, v0.getStart())
                                .attribute("SVTYPE", "INV")
                                .attribute("pctIdentity", wedit)
                                .attribute("prevBase", v0.getAttributeAsString("prevBase", "N"))
                                .attribute("nextBase", v1.getAttributeAsString("nextBase", "N"))
                                .rmAttribute("MATEID");

                        additions.add(vcb);

                        removals.add(vb0);
                        removals.add(vb1);
                    }
                }
            }
        }

        for (VariantContextBuilder addition : additions) {
            for (VariantContextBuilder call : calls) {
                if (addition.make().getStart() == call.make().getStart() || addition.make().getEnd() + 1 == call.make().getEnd()) {
                    removals.add(call);
                }
            }
        }

        callset.removeAll(removals);
        callset.addAll(additions);

        return callset;
    }

    private float pctSharedKmers(String a, String b) {
        Set<String> union = new HashSet<>(), ka = new HashSet<>(), kb = new HashSet<>();

        for (int i = 0; i <= a.length() - GRAPH.getKmerSize(); i++) {
            union.add(a.substring(i, i + GRAPH.getKmerSize()));
            ka.add(a.substring(i, i + GRAPH.getKmerSize()));
        }

        for (int i = 0; i <= b.length() - GRAPH.getKmerSize(); i++) {
            union.add(b.substring(i, i + GRAPH.getKmerSize()));
            kb.add(b.substring(i, i + GRAPH.getKmerSize()));
        }

        int numShared = 0;
        for (String kmer : union) {
            if (a.contains(kmer) && b.contains(kmer)) {
                numShared++;
            }
        }

        return 100.f * (float) numShared / (float) union.size();
    }

    private int numEdits(String[] a) {
        int numEdits = 0;
        for (int i = 0; i < a[0].length(); i++) {
            numEdits += a[0].charAt(i) == a[1].charAt(i) ? 0 : 1;
        }

        return numEdits;
    }

    private Set<VariantContextBuilder> mergeDoubleBreakpoints(String seq, Set<VariantContextBuilder> callset) {
        List<VariantContextBuilder> calls = new ArrayList<>(callset);

        if (calls.size() <= 1) {
            return callset;
        }

        List<VariantContextBuilder> bnds = new ArrayList<>();

        for (int i = 0; i < calls.size(); i++) {
            VariantContext vc = calls.get(i).make();
            if (vc.isSymbolicOrSV() && vc.getAttributeAsString("SVTYPE", "unknown").equals("BND")) {
                bnds.add(calls.get(i));
            }
        }

        Map<String, VariantContextBuilder> replacements = new HashMap<>();
        Set<Interval> removals = new HashSet<>();

        if (bnds.size() >= 4 && bnds.size() % 2 == 0) {
            for (int i = 0; i <= bnds.size() - 2; i += 2) {
                VariantContext outer0 = bnds.get(i).make();
                VariantContext inner0 = bnds.get(i+1).make();

                List<Triple<String, String, Pair<Integer, Integer>>> lps0 = (ArrayList<Triple<String, String, Pair<Integer, Integer>>>) outer0.getAttribute("lps");
                StringBuilder kmer0 = new StringBuilder();
                int pos0 = outer0.getAttributeAsInt("start", 0);
                while (kmer0.length() < GRAPH.getKmerSize() && pos0 >= 0) {
                    char c = getChildColumn(lps0, pos0);

                    if (c != '-' && c != ' ') {
                        kmer0.insert(0, c);
                    }
                }
                int q0 = getParentalRow(lps0, pos0);

                /*
                for (int q = 1; q < lps0.size(); q++) {
                    if (lps0.get(q).getLeft().equals(outer0.getAttributeAsString("targetName", ""))) {
                        int pos = outer0.getAttributeAsInt("start", 0);
                        while (kmer0.length() < GRAPH.getKmerSize() && pos >= 0) {
                            char c = lps0.get(q).getMiddle().charAt(pos);

                            if (c != '-') {
                                kmer0.insert(0, c);
                            }

                            pos--;
                        }

                        q0 = q;

                        break;
                    }
                }
                */

                for (int j = i + 2; j <= bnds.size() - 2; j += 2) {
                    VariantContext inner1 = bnds.get(j).make();
                    VariantContext outer1 = bnds.get(j+1).make();

                    List<Triple<String, String, Pair<Integer, Integer>>> lps1 = (ArrayList<Triple<String, String, Pair<Integer, Integer>>>) outer1.getAttribute("lps");
                    StringBuilder kmer1 = new StringBuilder();
                    int pos1 = outer1.getAttributeAsInt("start", 0);
                    while (kmer1.length() < GRAPH.getKmerSize() && pos1 < lps1.get(0).getMiddle().length()) {
                        char c = getChildColumn(lps1, pos1);

                        if (c != '-' && c != ' ') {
                            kmer1.append(c);
                        }
                    }
                    int q1 = getParentalRow(lps1, pos1);

                    /*
                    for (int q = 1; q < lps1.size(); q++) {
                        if (lps1.get(q).getLeft().equals(outer1.getAttributeAsString("targetName", ""))) {
                            int pos = outer1.getAttributeAsInt("start", 0);
                            while (kmer1.length() < GRAPH.getKmerSize() && pos < lps1.get(q).getMiddle().length()) {
                                char c = lps1.get(q).getMiddle().charAt(pos);

                                if (c != '-') {
                                    kmer1.append(c);
                                }

                                pos++;
                            }

                            q1 = q;

                            break;
                        }
                    }
                    */

                    String back0 = lps0.get(q0).getLeft().split(":")[0];
                    String back1 = lps1.get(q1).getLeft().split(":")[0];

                    if (back0.equals(back1)) {
                        for (Set<String> parentName : Arrays.asList(MOTHER, FATHER)) {
                            if (parentName.contains(back0)) {
                                TraversalEngine ef = new TraversalEngineFactory()
                                        .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                                        .traversalDirection(FORWARD)
                                        .combinationOperator(OR)
                                        .stoppingRule(DestinationStopper.class)
                                        .graph(GRAPH)
                                        .links(LINKS)
                                        .make();

                                TraversalEngine er = new TraversalEngineFactory()
                                        .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                                        .traversalDirection(REVERSE)
                                        .combinationOperator(OR)
                                        .stoppingRule(DestinationStopper.class)
                                        .graph(GRAPH)
                                        .links(LINKS)
                                        .make();

                                DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = ef.dfs(kmer0.toString(), kmer1.toString());
                                if (g == null || g.vertexSet().size() == 0) {
                                    g = er.dfs(kmer1.toString(), kmer0.toString());
                                }

                                if (g != null && g.vertexSet().size() > 0) {
                                    Set<String> parentalContigs = new HashSet<>();
                                    ConnectivityInspector<CortexVertex, CortexEdge> ci = new ConnectivityInspector<>(g);
                                    for (Set<CortexVertex> cvs : ci.connectedSets()) {
                                        List<CortexVertex> w = new ArrayList<>();

                                        for (CortexVertex cv : cvs) {
                                            List<CortexVertex> wa = TraversalUtils.toWalk(g, cv.getKmerAsString(), g.edgeSet().iterator().next().getColor());

                                            if (wa.size() == w.size()) { break; }
                                            else if (wa.size() > w.size()) {
                                                w = wa;
                                            }
                                        }

                                        if (w.size() > 0) {
                                            parentalContigs.add(TraversalUtils.toContig(w));
                                        }
                                    }

                                    String childContig = seq.substring(seq.indexOf(kmer0.toString()) + GRAPH.getKmerSize(), seq.indexOf(kmer1.toString(), seq.indexOf(kmer0.toString())));

                                    for (String parentalContig : parentalContigs) {
                                        String inverted = SequenceUtils.reverseComplement(parentalContig.substring(GRAPH.getKmerSize(), parentalContig.length() - GRAPH.getKmerSize()));

                                        SmithWaterman sw = new SmithWaterman();
                                        String[] a = sw.getAlignment(childContig, inverted);

                                        int edits = 0;
                                        for (int l = 0; l < a[0].length(); l++) {
                                            if (a[0].charAt(l) != a[1].charAt(l)) {
                                                edits++;
                                            }
                                        }

                                        float pctIdentity = 100.0f * (((float) (a[0].length() - edits)) / ((float) a[0].length()));

                                        if (pctIdentity >= 90.0) {
                                            List<Allele> alleles = Arrays.asList(Allele.create(parentalContig, true), Allele.create(childContig));

                                            VariantContextBuilder vcb = new VariantContextBuilder(outer0)
                                                    .alleles(alleles)
                                                    .computeEndFromAlleles(alleles, outer0.getStart())
                                                    .attribute("SVTYPE", "INV")
                                                    .attribute("pctIdentity", pctIdentity)
                                                    .attribute("prevBase", outer0.getAttributeAsString("prevBase", "N"))
                                                    .attribute("nextBase", outer1.getAttributeAsString("nextBase", "N"))
                                                    .rmAttribute("MATEID");

                                            replacements.put(outer0.getID(), vcb);
                                            replacements.put(inner0.getID(), null);
                                            replacements.put(inner1.getID(), null);
                                            replacements.put(outer1.getID(), null);

                                            removals.add(new Interval(outer0.getContig(), outer0.getStart(), outer0.getStart()));
                                            removals.add(new Interval(inner0.getContig(), inner0.getStart(), inner0.getStart()));
                                            removals.add(new Interval(inner1.getContig(), inner1.getStart(), inner1.getStart()));
                                            removals.add(new Interval(outer1.getContig(), outer1.getStart(), outer1.getStart()));
                                        }
                                    }
                                } else {

                                }
                            }
                        }
                    }
                }
            }
        }

        Set<VariantContextBuilder> newcalls = new LinkedHashSet<>();
        for (VariantContextBuilder vcb : calls) {
            if (!vcb.make().isSymbolic() && removals.contains(new Interval(vcb.make().getContig(), vcb.make().getStart(), vcb.make().getStart()))) {
                // ignore
            } else {
                if (!replacements.containsKey(vcb.make().getID())) {
                    newcalls.add(vcb);
                } else {
                    if (replacements.get(vcb.make().getID()) != null) {
                        newcalls.add(replacements.get(vcb.make().getID()));
                    }
                }
            }
        }

        return newcalls;
    }

    private List<VariantContextBuilder> callLargeBubbles(List<Triple<String, String, Pair<Integer, Integer>>> lps, List<Pair<Integer, Integer>> nrs, Map<String, String> targets, String contigName, int sectionStart, int sectionStop) {
        List<VariantContextBuilder> vcbs = new ArrayList<>();

        for (Pair<Integer, Integer> nr : nrs) {
            for (int i = nr.getFirst(); i <= nr.getSecond(); i++) {
                if (isRecomb(lps, i)) {
                    Pair<Integer, Integer> partners = recombPartners(lps, i);

                    String name0 = lps.get(partners.getFirst()).getLeft();
                    String name1 = lps.get(partners.getSecond()).getLeft();

                    if (name0.equals(name1)) {
                        String target = targets.get(lps.get(partners.getFirst()).getLeft());

                        int start = lps.get(partners.getFirst()).getRight().getSecond() + 1;
                        int stop = lps.get(partners.getSecond()).getRight().getFirst();

                        if (stop > start) {
                            int variantStart = sectionStart + i;
                            int variantStop = sectionStart + i + 1;

                            char prevBase = Character.toUpperCase(getParentalColumn(lps, i));
                            char nextBase = Character.toUpperCase(getParentalColumn(lps, i + 1));

                            String subtarget = target.substring(start, stop);

                            List<Allele> alleles = Arrays.asList(Allele.create(String.valueOf(prevBase), true), Allele.create(String.valueOf(prevBase) + subtarget));

                            String varBackground = lps.get(partners.getFirst()).getLeft().split(":")[0];

                            int childLeft;
                            int numLeft = 0;
                            for (childLeft = nr.getFirst(); childLeft > 0 && numLeft <= GRAPH.getKmerSize(); childLeft--) {
                                if (getChildColumn(lps, childLeft) != '-') {
                                    numLeft++;
                                }
                            }

                            int childRight;
                            int numRight = 0;
                            for (childRight = nr.getSecond(); childRight < lps.get(0).getMiddle().length() && numRight <= GRAPH.getKmerSize(); childRight++) {
                                if (getChildColumn(lps, childRight) != '-') {
                                    numRight++;
                                }
                            }

                            String childHap = lps.get(0).getMiddle().substring(childLeft, childRight).replaceAll("-", "");

                            VariantContextBuilder vcb = new VariantContextBuilder()
                                    .noID()
                                    .noGenotypes()
                                    .chr(contigName)
                                    .alleles(alleles)
                                    .start(variantStart)
                                    .computeEndFromAlleles(alleles, sectionStart + i)
                                    .attribute("start", i)
                                    .attribute("stop", i + 1)
                                    .attribute("sectionStart", sectionStart)
                                    .attribute("sectionStop", sectionStop)
                                    .attribute("variantStart", variantStart)
                                    .attribute("variantStop", variantStop)
                                    .attribute("prevBase", prevBase)
                                    .attribute("nextBase", nextBase)
                                    .attribute("CHILD_HAP", childHap)
                                    .attribute("PARTITION_NAME", contigName)
                                    .attribute("BACKGROUND", varBackground);

                            vcbs.add(vcb);
                        }
                    }
                }
            }
        }

        return vcbs;
    }

    private List<VariantContextBuilder> callBreakpoints(List<Triple<String, String, Pair<Integer, Integer>>> lps, List<Pair<Integer, Integer>> nrs, String contigName, int sectionStart, int sectionStop) {
        List<VariantContextBuilder> vcbs = new ArrayList<>();

        for (Pair<Integer, Integer> nr : nrs) {
            for (int i = nr.getFirst(); i <= nr.getSecond(); i++) {
                if (isRecomb(lps, i)) {
                    Pair<Integer, Integer> partners = recombPartners(lps, i);

                    String name0 = lps.get(partners.getFirst()).getLeft();
                    String name1 = lps.get(partners.getSecond()).getLeft();

                    if (!name0.equals(name1)) {
                        int prevPos = i, nextPos = i + 1;

                        StringBuilder nextIns = new StringBuilder();
                        while (getParentalColumn(lps, prevPos) == '-') {
                            nextIns.insert(0, getChildColumn(lps, prevPos));
                            prevPos--;
                        }
                        nextIns.insert(0, getChildColumn(lps, prevPos));
                        char prevBase = getChildColumn(lps, prevPos);

                        StringBuilder prevIns = new StringBuilder();
                        while (getParentalColumn(lps, nextPos) == '-') {
                            prevIns.append(getChildColumn(lps, nextPos));
                            nextPos++;
                        }
                        prevIns.append(getChildColumn(lps, nextPos));
                        char nextBase = getChildColumn(lps, nextPos);

                        List<Allele> a0 = Arrays.asList(Allele.create(String.valueOf(prevBase), true), Allele.create("]" + name1 + ":" + nextPos + "]" + nextIns.toString()));
                        List<Allele> a1 = Arrays.asList(Allele.create(String.valueOf(nextBase), true), Allele.create(prevIns.toString() + "[" + name0 + ":" + prevPos + "["));

                        String mate0 = String.format("bnd_%s_%d", contigName, sectionStart + prevPos);
                        String mate1 = String.format("bnd_%s_%d", contigName, sectionStart + nextPos);

                        String varBackground0 = lps.get(partners.getFirst()).getLeft().split(":")[0];
                        String varBackground1 = lps.get(partners.getSecond()).getLeft().split(":")[0];

                        int childLeft;
                        int numLeft = 0;
                        for (childLeft = nr.getFirst(); childLeft > 0 && numLeft <= GRAPH.getKmerSize(); childLeft--) {
                            if (getChildColumn(lps, childLeft) != '-') {
                                numLeft++;
                            }
                        }

                        int childRight;
                        int numRight = 0;
                        for (childRight = nr.getSecond(); childRight < lps.get(0).getMiddle().length() && numRight <= GRAPH.getKmerSize(); childRight++) {
                            if (getChildColumn(lps, childRight) != '-') {
                                numRight++;
                            }
                        }

                        String childHap = lps.get(0).getMiddle().substring(childLeft, childRight).replaceAll("-", "");

                        VariantContextBuilder vcb0 = new VariantContextBuilder()
                                .id(mate0)
                                .noGenotypes()
                                .chr(contigName)
                                .alleles(a0)
                                .start(sectionStart + prevPos)
                                .stop(sectionStart + prevPos)
                                .attribute("start", prevPos)
                                .attribute("stop", prevPos + 1)
                                .attribute("sectionStart", sectionStart)
                                .attribute("sectionStop", sectionStop)
                                .attribute("variantStart", sectionStart + prevPos)
                                .attribute("variantStop", sectionStart + prevPos)
                                .attribute("prevBase", prevBase)
                                .attribute("nextBase", nextBase)
                                .attribute("targetName", lps.get(partners.getFirst()).getLeft())
                                .attribute("targetStart", lps.get(partners.getFirst()).getRight().getFirst())
                                .attribute("targetStop", lps.get(partners.getFirst()).getRight().getSecond())
                                .attribute("CHILD_HAP", childHap)
                                .attribute("PARTITION_NAME", contigName)
                                .attribute("BACKGROUND", varBackground0)
                                .attribute("SVTYPE", "BND")
                                .attribute("MATEID", mate1);

                        VariantContextBuilder vcb1 = new VariantContextBuilder()
                                .id(mate1)
                                .noGenotypes()
                                .chr(contigName)
                                .alleles(a1)
                                .start(sectionStart + nextPos)
                                .stop(sectionStart + nextPos)
                                .attribute("start", nextPos)
                                .attribute("stop", nextPos + 1)
                                .attribute("sectionStart", sectionStart)
                                .attribute("sectionStop", sectionStop)
                                .attribute("variantStart", sectionStart + nextPos)
                                .attribute("variantStop", sectionStart + nextPos)
                                .attribute("prevBase", prevBase)
                                .attribute("nextBase", nextBase)
                                .attribute("targetName", lps.get(partners.getSecond()).getLeft())
                                .attribute("targetStart", lps.get(partners.getSecond()).getRight().getFirst())
                                .attribute("targetStop", lps.get(partners.getSecond()).getRight().getSecond())
                                .attribute("CHILD_HAP", childHap)
                                .attribute("PARTITION_NAME", contigName)
                                .attribute("BACKGROUND", varBackground1)
                                .attribute("SVTYPE", "BND")
                                .attribute("MATEID", mate0);

                        vcbs.add(vcb0);
                        vcbs.add(vcb1);
                    }
                }
            }
        }

        return vcbs;
    }

    private List<VariantContextBuilder> callSmallBubbles(List<Triple<String, String, Pair<Integer, Integer>>> lps, List<Pair<Integer, Integer>> nrs, String contigName, int sectionStart, int sectionStop) {
        List<VariantContextBuilder> vcbs = new ArrayList<>();

        for (Pair<Integer, Integer> nr : nrs) {
            VariantContextBuilder vcb = new VariantContextBuilder();
            int start = nr.getFirst() - 1;
            char prevBase = getChildColumn(lps, start);
            int prevRow = getParentalRow(lps, start);
            StringBuilder cBuilder = null;
            StringBuilder pBuilder = null;

            for (int i = nr.getFirst(); i <= nr.getSecond(); i++) {
                if (Character.toUpperCase(getChildColumn(lps, i)) == Character.toUpperCase(getParentalColumn(lps, i)) || i == getNumColumns(lps) - 1) {
                    if (cBuilder != null) {
                        if (i == getNumColumns(lps) - 1) {
                            if (getChildColumn(lps, i) != '-') { cBuilder.append(getChildColumn(lps, i)); }
                            if (getParentalColumn(lps, i) != '-') { pBuilder.append(getParentalColumn(lps, i)); }

                            cBuilder.append(".");
                        }

                        boolean isSymbolicStart = cBuilder.length() > 0 && cBuilder.charAt(0) == '.';
                        boolean isSymbolicEnd = cBuilder.length() > 0 && cBuilder.charAt(cBuilder.length() - 1) == '.';

                        int variantStart = sectionStart + start;
                        int variantStop = sectionStart + i;
                        char nextBase = i == getNumColumns(lps) - 1 ? 'N' : getChildColumn(lps, i);
                        int nextRow = getParentalRow(lps, i);

                        if (cBuilder.length() == pBuilder.length() && cBuilder.length() == 1) {
                            variantStart++;
                            variantStop--;
                        } else {
                            if (!isSymbolicStart) {
                                cBuilder.insert(0, prevBase);
                                pBuilder.insert(0, prevBase);
                            } else {
                                variantStart = variantStop;
                                start = i;
                                cBuilder.append(nextBase);
                                pBuilder.append(nextBase);
                            }
                        }

                        int childLeft;
                        int numLeft = 0;
                        for (childLeft = nr.getFirst(); childLeft > 0 && numLeft <= GRAPH.getKmerSize(); childLeft--) {
                            if (getChildColumn(lps, childLeft) != '-') {
                                numLeft++;
                            }
                        }

                        int childRight;
                        int numRight = 0;
                        for (childRight = nr.getSecond(); childRight < lps.get(0).getMiddle().length() && numRight <= GRAPH.getKmerSize(); childRight++) {
                            if (getChildColumn(lps, childRight) != '-') {
                                numRight++;
                            }
                        }

                        String childHap = lps.get(0).getMiddle().substring(childLeft, childRight).replaceAll("-", "");

                        int row = prevRow == 0 ? nextRow : prevRow;

                        List<Allele> alleleStrings = Arrays.asList(Allele.create(pBuilder.toString(), true), Allele.create(cBuilder.toString()));

                        String varBackground = row > 0 ? lps.get(row).getLeft().split(":")[0] : "unknown";

                        vcb
                                .noID()
                                .noGenotypes()
                                .chr(contigName)
                                .alleles(alleleStrings)
                                .start(variantStart)
                                .attribute("start", start)
                                .attribute("stop", i)
                                .attribute("sectionStart", sectionStart)
                                .attribute("sectionStop", sectionStop)
                                .attribute("variantStart", variantStart)
                                .attribute("variantStop", variantStop)
                                .attribute("prevBase", prevBase)
                                .attribute("nextBase", nextBase)
                                .attribute("CHILD_HAP", childHap)
                                .attribute("PARTITION_NAME", contigName)
                                .attribute("BACKGROUND", varBackground)
                        ;

                        if (isSymbolicStart || isSymbolicEnd) {
                            vcb.stop(variantStop)
                               .attribute("SVTYPE", "BND");

                        } else {
                            vcb.computeEndFromAlleles(alleleStrings, variantStart);
                        }

                        vcbs.add(vcb);
                    }

                    vcb = new VariantContextBuilder();
                    prevBase = getChildColumn(lps, i);
                    start = i;
                    cBuilder = null;
                    pBuilder = null;
                } else {
                    if (cBuilder == null) { cBuilder = new StringBuilder(); }
                    if (pBuilder == null) { pBuilder = new StringBuilder(); }

                    if (i == 0) { cBuilder.insert(0, "."); }

                    if (getChildColumn(lps, i) != '-') {
                        cBuilder.append(getChildColumn(lps, i));
                    }
                    if (getParentalColumn(lps, i) != '-') {
                        pBuilder.append(getParentalColumn(lps, i));
                    }
                }
            }
        }

        return vcbs;
    }

    private void writeVariants(Set<CanonicalKmer> rois, Set<VariantContext> svcs, VariantContextWriter vcw) {
        Map<CanonicalKmer, String> acct = new TreeMap<>();
        for (CanonicalKmer ck : rois) {
            acct.put(ck, "absent");
        }

        int variantId = 0;
        for (VariantContext vc : svcs) {
            String id = String.format("CC%d", variantId);

            vcw.add(new VariantContextBuilder(vc)
                .rmAttribute("novels")
                .attribute("CALL_ID", variantId)
                //.id(id)
                .make()
            );

            for (String sk : vc.getAttributeAsString("novels", "").split(",")) {
                if (sk.length() > 0) {
                    CanonicalKmer ck = new CanonicalKmer(sk);

                    if (acct.containsKey(ck)) {
                        acct.put(ck, id);
                    }
                }
            }

            variantId++;
        }

        vcw.close();

        for (CanonicalKmer ck : acct.keySet()) {
            aout.println(Joiner.on("\t").join(ck, acct.get(ck)));
        }
    }

    @NotNull
    private VariantContextWriter buildVariantWriter(SAMSequenceDictionary ssd) {
        VariantContextWriter vcw = new VariantContextWriterBuilder()
                .setOutputFile(out)
                .setOutputFileType(VCF)
                .setOption(Options.DO_NOT_WRITE_GENOTYPES)
                .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                .unsetOption(Options.INDEX_ON_THE_FLY)
                .build();

        VCFHeader vcfHeader = new VCFHeader();
        vcfHeader.setSequenceDictionary(ssd);
        vcw.writeHeader(vcfHeader);
        return vcw;
    }

    @NotNull
    private Set<VariantContext> buildVariantSorter(SAMSequenceDictionary ssd) {
        Map<String, Integer> sid = new HashMap<>();
        for (int i = 0; i < ssd.getSequences().size(); i++) {
            sid.put(ssd.getSequence(i).getSequenceName(), i);
        }

        return new TreeSet<>((v1, v2) -> {
            if (v1 != null && v2 != null) {
                int sid0 = sid.getOrDefault(v1.getContig(), 0);
                int sid1 = sid.getOrDefault(v2.getContig(), 0);

                if (sid0 != sid1) { return sid0 < sid1 ? -1 : 1; }
                if (v1.getStart() != v2.getStart()) { return v1.getStart() < v2.getStart() ? -1 : 1; }
                if (v1.isSymbolic() != v2.isSymbolic()) { return v1.isSymbolic() ? 1 : -1; }
            }

            return 0;
        });
    }

    @NotNull
    private Set<VariantContextBuilder> buildVariantContextBuilderSorter(SAMSequenceDictionary ssd) {
        Map<String, Integer> sid = new HashMap<>();
        for (int i = 0; i < ssd.getSequences().size(); i++) {
            sid.put(ssd.getSequence(i).getSequenceName(), i);
        }

        return new TreeSet<>((vb1, vb2) -> {
            VariantContext v1 = vb1.make();
            VariantContext v2 = vb2.make();

            int sid0 = sid.getOrDefault(v1.getContig(), 0);
            int sid1 = sid.getOrDefault(v2.getContig(), 0);

            if (sid0 != sid1) { return sid0 < sid1 ? -1 : 1; }
            if (v1.getStart() != v2.getStart()) { return v1.getStart() < v2.getStart() ? -1 : 1; }
            if (v1.isSymbolic() != v2.isSymbolic()) { return v1.isSymbolic() ? 1 : -1; }

            return 0;
        });
    }

    @NotNull
    private SAMSequenceDictionary buildMergedSequenceDictionary(List<ReferenceSequence> rseqs) {
        List<SAMSequenceRecord> ssrs = new ArrayList<>();
        for (String id : BACKGROUNDS.keySet()) {
            IndexedReference ir = BACKGROUNDS.get(id);
            ssrs.addAll(ir.getReferenceSequence().getSequenceDictionary().getSequences());
            ssrs.add(new SAMSequenceRecord(id + "_unknown", rseqs.size()));
        }
        return new SAMSequenceDictionary(ssrs);
    }

    @NotNull
    private List<ReferenceSequence> loadPartitions() {
        List<ReferenceSequence> rseqs = new ArrayList<>();
        ReferenceSequence rseq;
        while ((rseq = PARTITIONS.nextSequence()) != null) {
            String[] name = rseq.getName().split(" ");
            if (PARTITION_NAMES == null || PARTITION_NAMES.contains(name[0])) {
                rseqs.add(rseq);
            }
        }
        return rseqs;
    }

    private List<SAMRecord> sortAlignments(String background, String target) {
        if (!BACKGROUNDS.containsKey(background)) {
            return new ArrayList<>();
        }

        List<SAMRecord> a = BACKGROUNDS.get(background).align(target);
        //a.removeIf(s -> s.getMappingQuality() == 0);

        a.sort((s1, s2) -> {
            int s1length = s1.getAlignmentEnd() - s1.getAlignmentStart();
            int nm1 = s1.getIntegerAttribute("NM");
            int mq1 = s1.getMappingQuality();

            int s2length = s2.getAlignmentEnd() - s2.getAlignmentStart();
            int nm2 = s2.getIntegerAttribute("NM");
            int mq2 = s1.getMappingQuality();

            if (s1length != s2length) {
                return s1length > s2length ? -1 : 1;
            }

            if (nm1 != nm2) {
                return nm1 < nm2 ? -1 : 1;
            }

            if (mq1 != mq2) {
                return mq1 > mq2 ? -1 : 1;
            }

            return 0;
        });

        return a;
    }

    private Triple<Integer, Integer, String> trimQuery(List<CortexVertex> ws, Map<String, String> targets, Set<CanonicalKmer> rois) {
        int firstIndex = Integer.MAX_VALUE, lastIndex = 0;
        int firstNovel = -1, lastNovel = -1;

        Map<CanonicalKmer, Integer> pos = new HashMap<>();
        for (int i = 0; i < ws.size(); i++) {
            pos.put(ws.get(i).getCanonicalKmer(), i);

            if (rois.contains(ws.get(i).getCanonicalKmer())) {
                if (firstNovel == -1) {
                    firstNovel = i;
                }

                lastNovel = i;
            }
        }

        for (String target : targets.values()) {
            for (int i = 0; i <= target.length() - GRAPH.getKmerSize(); i++) {
                String sk = target.substring(i, i + GRAPH.getKmerSize());
                CanonicalKmer ck = new CanonicalKmer(sk);

                if (pos.containsKey(ck)) {
                    int index = pos.get(ck);

                    if (index < firstIndex) { firstIndex = index; }
                    if (index > lastIndex) { lastIndex = index; }
                }
            }
        }

        if (firstNovel < firstIndex) { firstIndex = firstNovel; }
        if (lastNovel > lastIndex) { lastIndex = lastNovel; }

        return Triple.of(firstIndex, lastIndex, TraversalUtils.toContig(ws.subList(firstIndex, lastIndex)));
    }

    private int getNumColumns(List<Triple<String, String, Pair<Integer, Integer>>> lps) {
        return lps.get(0).getMiddle().length();
    }

    private Pair<Integer, Integer> recombPartners(List<Triple<String, String, Pair<Integer, Integer>>> lps, int column) {
        if (lps.size() > 2) {
            for (int i = 1; i < lps.size(); i++) {
                if (column == lps.get(i).getMiddle().length() - 1 && lps.get(i).getMiddle().charAt(column) != ' ' &&
                    column + 1 < lps.get(i + 1).getMiddle().length() && lps.get(i + 1).getMiddle().charAt(column) == ' ' && lps.get(i + 1).getMiddle().charAt(column + 1) != ' ') {
                    return Pair.create(i, i+1);
                }
            }
        }

        return Pair.create(-1, -1);
    }

    private boolean isRecomb(List<Triple<String, String, Pair<Integer, Integer>>> lps, int column) {
        if (lps.size() > 2) {
            for (int i = 1; i < lps.size() - 1; i++) {
                if (column == lps.get(i).getMiddle().length() - 1 && lps.get(i).getMiddle().charAt(column) != ' ' &&
                    column + 1 < lps.get(i + 1).getMiddle().length() && lps.get(i + 1).getMiddle().charAt(column) == ' ' && lps.get(i + 1).getMiddle().charAt(column + 1) != ' ') {

                    return true;
                }
            }
        }

        return false;
    }

    private char getChildColumn(List<Triple<String, String, Pair<Integer, Integer>>> lps, int column) {
        if (column >= 0 && column < getNumColumns(lps)) {
            if (column < lps.get(0).getMiddle().length() && lps.get(0).getMiddle().charAt(column) != ' ') {
                return lps.get(0).getMiddle().charAt(column);
            }
        }

        return 'N';
    }

    private char getParentalColumn(List<Triple<String, String, Pair<Integer, Integer>>> lps, int column) {
        if (column >= 0 && column < getNumColumns(lps)) {
            for (int i = 1; i < lps.size(); i++) {
                if (column < lps.get(i).getMiddle().length() && lps.get(i).getMiddle().charAt(column) != ' ') {
                    return lps.get(i).getMiddle().charAt(column);
                }
            }
        }

        return 'N';
    }

    private int getParentalRow(List<Triple<String, String, Pair<Integer, Integer>>> lps, int column) {
        if (column >= 0 && column < getNumColumns(lps)) {
            for (int i = 1; i < lps.size(); i++) {
                if (column < lps.get(i).getMiddle().length() && lps.get(i).getMiddle().charAt(column) != ' ') {
                    return i;
                }
            }
        }

        return 0;
    }

    private String makeNoveltyTrack(Set<CanonicalKmer> rois, List<Triple<String, String, Pair<Integer, Integer>>> lps, boolean expand) {
        int maxLength = 0;
        for (Triple<String, String, Pair<Integer, Integer>> lp : lps) {
            String name = String.format("%s (%d-%d)", lp.getLeft(), lp.getRight().getFirst(), lp.getRight().getSecond());
            maxLength = Math.max(maxLength, name.length());
        }

        String query = lps.get(0).getMiddle().replaceAll("[- ]", "");
        StringBuilder sb = new StringBuilder(StringUtil.repeatCharNTimes(' ', query.length() + 1));
        for (int i = 0; i <= query.length() - GRAPH.getKmerSize(); i++) {
            CanonicalKmer ck = new CanonicalKmer(query.substring(i, i + GRAPH.getKmerSize()));

            if (rois.contains(ck)) {
                for (int j = i; j < i + GRAPH.getKmerSize(); j++) {
                    sb.setCharAt(j, '*');
                }
            }
        }

        for (int i = 0; i < getNumColumns(lps); i++) {
            if (getChildColumn(lps, i) == '-') {
                sb.insert(i, sb.charAt(i) == '*' ? '*' : ' ');
            }
        }

        if (expand) {
            for (int i = 1; i < getNumColumns(lps); i++) {
                if (sb.charAt(i) == '*') {
                    if (sb.charAt(i - 1) != '*' && getParentalColumn(lps, i - 1) == '-') {
                        for (int j = i - 1; j >= 0 && getParentalColumn(lps, j) == '-'; j--) {
                            sb.setCharAt(j, '*');
                        }
                    }

                    if (sb.charAt(i + 1) != '*' && getParentalColumn(lps, i + 1) == '-') {
                        for (int j = i + 1; j < getNumColumns(lps) && getParentalColumn(lps, j) == '-'; j++) {
                            sb.setCharAt(j, '*');
                        }
                    }
                }
            }
        }

        return String.format("%" + maxLength + "s %s", "novel", sb.toString());
    }

    private List<Pair<Integer, Integer>> getNoveltyRegions(Set<CanonicalKmer> rois, List<Triple<String, String, Pair<Integer, Integer>>> lps, boolean expand) {
        String noveltyTrack = makeNoveltyTrack(rois, lps, expand);
        noveltyTrack = noveltyTrack.replaceAll("^\\s+novel ", "");

        List<Pair<Integer, Integer>> regions = new ArrayList<>();
        int start = -1;
        int stop = noveltyTrack.length() - 1;
        for (int i = 0; i < noveltyTrack.length(); i++) {
            if (noveltyTrack.charAt(i) == '*') {
                if (start == -1) { start = i; }
                stop = i;
            } else {
                if (start >= 0) {
                    regions.add(Pair.create(start, stop));
                    start = -1;
                    stop = noveltyTrack.length() - 1;
                }
            }
        }

        if (start >= 0) {
            regions.add(Pair.create(start, stop));
        }

        return regions;
    }

    private Map<String, String> fasterAssembleCandidateHaplotypes(List<CortexVertex> ws, Set<String> parentName) {
        List<Integer> colors = GRAPH.getColorsForSampleNames(parentName);

        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = new DirectedWeightedPseudograph<>(CortexEdge.class);

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColors(colors)
                .traversalDirection(BOTH)
                .combinationOperator(OR)
                .stoppingRule(ContigStopper.class)
                .maxBranchLength(ws.size())
                .graph(GRAPH)
                .links(LINKS)
                .make();

        for (int i = 0; i < ws.size(); i++) {
            boolean hasCoverage = false;
            for (int c : colors) {
                hasCoverage |= ws.get(i).getCortexRecord().getCoverage(c) > 0;
            }

            if (hasCoverage && TraversalUtils.findVertex(g, ws.get(i).getKmerAsString()) == null) {
                DirectedWeightedPseudograph<CortexVertex, CortexEdge> gs = e.dfs(ws.get(i).getKmerAsString());

                if (gs != null && gs.vertexSet().size() > 0) {
                    Graphs.addGraph(g, gs);
                }
            }
        }

        Set<CortexVertex> inEnds = getCloseableGraphEnds(colors, g, false);
        Set<CortexVertex> outEnds = getCloseableGraphEnds(colors, g, true);
        closeGaps(colors, g, inEnds, outEnds);
        extendFlanks(colors, g, inEnds, outEnds);

        Map<String, String> targets = new HashMap<>();
        if (g.edgeSet().size() > 0) {
            int repColor = g.edgeSet().iterator().next().getColor();

            ConnectivityInspector<CortexVertex, CortexEdge> ci = new ConnectivityInspector<>(g);
            List<List<CortexVertex>> walks = new ArrayList<>();

            for (Set<CortexVertex> cs : ci.connectedSets()) {
                List<CortexVertex> w = new ArrayList<>();

                for (CortexVertex cv : cs) {
                    List<CortexVertex> wa = TraversalUtils.toWalk(g, cv.getKmerAsString(), repColor);

                    if (wa.size() == w.size()) {
                        break;
                    } else if (wa.size() > w.size()) {
                        w = wa;
                    }
                }

                if (w.size() > 0) {
                    walks.add(w);
                }
            }

            Set<CanonicalKmer> indices = new HashSet<>();
            for (CortexVertex cv : ws) {
                indices.add(cv.getCanonicalKmer());
            }

            Set<String> contigs = new HashSet<>();
            for (List<CortexVertex> w : walks) {
                int actualStart = Integer.MAX_VALUE, actualEnd = -1;
                int shared = 0;
                for (int i = 0; i < w.size(); i++) {
                    if (indices.contains(w.get(i).getCanonicalKmer())) {
                        shared++;

                        if (i < actualStart) { actualStart = i; }
                        if (i > actualEnd) { actualEnd = i; }
                    }
                }

                if (actualStart == Integer.MAX_VALUE) { actualStart = 0; }
                if (actualEnd == -1 || actualEnd == actualStart) { actualEnd = w.size() - 1; }

                if (shared > 0) {
                    String contig = TraversalUtils.toContig(w.subList(actualStart, actualEnd));

                    contigs.add(contig);
                }
            }

            int i = 0;
            for (String contig : contigs) {
                String id = String.format("%s:%s_unknown:%s_contig%d_%s", parentName.iterator().next(), parentName.iterator().next(), parentName.iterator().next(), i, "fastasm");

                if (contig.length() > 0) {
                    targets.put(id, contig);

                    i++;
                }
            }
        }

        return targets;
    }

    private void closeGaps(List<Integer> colors, DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, Set<CortexVertex> inEnds, Set<CortexVertex> outEnds) {
        TraversalEngine ef = new TraversalEngineFactory()
                .traversalColors(colors)
                .traversalDirection(FORWARD)
                .combinationOperator(OR)
                .stoppingRule(DestinationStopper.class)
                .graph(GRAPH)
                .links(LINKS)
                .make();

        TraversalEngine er = new TraversalEngineFactory()
                .traversalColors(colors)
                .traversalDirection(REVERSE)
                .combinationOperator(OR)
                .stoppingRule(DestinationStopper.class)
                .graph(GRAPH)
                .links(LINKS)
                .make();

        for (CortexVertex ie : inEnds) {
            for (CortexVertex oe : outEnds) {
                DirectedWeightedPseudograph<CortexVertex, CortexEdge> gg = ef.dfs(ie.getKmerAsString(), oe.getKmerAsString());
                if (gg == null || gg.vertexSet().size() == 0) {
                    gg = er.dfs(oe.getKmerAsString(), ie.getKmerAsString());
                }

                if (gg != null && gg.vertexSet().size() > 0) {
                    Graphs.addGraph(g, gg);
                }
            }
        }
    }

    private void extendFlanks(List<Integer> colors, DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, Set<CortexVertex> inEnds, Set<CortexVertex> outEnds) {
        TraversalEngine eb = new TraversalEngineFactory()
                .traversalColors(colors)
                .traversalDirection(BOTH)
                .combinationOperator(OR)
                .stoppingRule(ContigStopper.class)
                .maxBranchLength(500)
                .graph(GRAPH)
                .links(LINKS)
                .make();

        for (Set<CortexVertex> cvs : Arrays.asList(inEnds, outEnds)) {
            for (CortexVertex cv : cvs) {
                DirectedWeightedPseudograph<CortexVertex, CortexEdge> gg = eb.dfs(cv.getKmerAsString());

                if (gg != null && gg.vertexSet().size() > 0) {
                    Graphs.addGraph(g, gg);
                }
            }
        }
    }

    @NotNull
    private Set<CortexVertex> getCloseableGraphEnds(List<Integer> colors, DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, boolean outgoing) {
        Set<CortexVertex> ends = new HashSet<>();
        if (g.edgeSet().size() > 0) {
            for (CortexVertex cv : g.vertexSet()) {
                if (outgoing) {
                    if (g.outDegreeOf(cv) == 0) {
                        ends.addAll(Graphs.predecessorListOf(g, cv));
                    }
                } else {
                    if (g.inDegreeOf(cv) == 0) {
                        ends.addAll(Graphs.successorListOf(g, cv));
                    }
                }
            }
        }

        TraversalEngine ef = new TraversalEngineFactory()
                .traversalColors(colors)
                .traversalDirection(FORWARD)
                .combinationOperator(OR)
                .stoppingRule(ContigStopper.class)
                .maxBranchLength(10)
                .graph(GRAPH)
                .links(LINKS)
                .make();

        TraversalEngine er = new TraversalEngineFactory()
                .traversalColors(colors)
                .traversalDirection(REVERSE)
                .combinationOperator(OR)
                .stoppingRule(ContigStopper.class)
                .maxBranchLength(10)
                .graph(GRAPH)
                .links(LINKS)
                .make();

        Set<CortexVertex> endsToRemove = new HashSet<>();
        for (CortexVertex e0 : ends) {
            String fw = e0.getKmerAsString();

            for (CortexVertex e1 : ends) {
                if (!e0.equals(e1) && !endsToRemove.contains(e0) && !endsToRemove.contains(e1)) {
                    String rc = SequenceUtils.reverseComplement(e1.getKmerAsString());

                    DirectedWeightedPseudograph<CortexVertex, CortexEdge> gf = ef.dfs(fw, rc);
                    DirectedWeightedPseudograph<CortexVertex, CortexEdge> gr = er.dfs(rc, fw);

                    if ((gf != null && gf.vertexSet().size() > 0) || (gr != null && gr.vertexSet().size() > 0)) {
                        endsToRemove.add(e0);
                        endsToRemove.add(e1);
                    }
                }
            }
        }

        ends.removeAll(endsToRemove);

        return ends;
    }

    private Set<CanonicalKmer> loadRois(CortexGraph rg) {
        Set<CanonicalKmer> rois = new HashSet<>();

        for (CortexRecord rr : rg) {
            rois.add(rr.getCanonicalKmer());
        }

        return rois;
    }

    private List<CortexVertex> loadChildWalk(ReferenceSequence seq, CortexGraph graph) {
        List<CortexVertex> w = new ArrayList<>();

        Map<String, Integer> seenCount = new HashMap<>();

        String contig = seq.getBaseString();
        for (int i = 0; i <= contig.length() - graph.getKmerSize(); i++) {
            String sk = contig.substring(i, i + graph.getKmerSize());

            if (!seenCount.containsKey(sk)) {
                seenCount.put(sk, 0);
            } else {
                ContainerUtils.increment(seenCount, sk);
            }

            w.add(new CortexVertexFactory()
                    .bases(sk)
                    .record(graph.findRecord(sk))
                    .copyIndex(seenCount.get(sk))
                    .make());
        }

        return w;
    }

    private List<Triple<Integer, Integer, List<CortexVertex>>> sectionContig(Set<CanonicalKmer> rois, List<CortexVertex> w, int window, int novelDistanceSplit) {
        List<Pair<Integer, Integer>> regions = getRegions(rois, w);

        if (regions.size() > 0) {
            int subcontigStart = regions.get(0).getFirst() - window;
            if (subcontigStart < 0) {
                subcontigStart = 0;
            }

            int subcontigStop = regions.get(regions.size() - 1).getSecond() + window;
            if (subcontigStop >= w.size()) {
                subcontigStop = w.size() - 1;
            }

            List<Pair<Integer, Integer>> sections = new ArrayList<>();
            for (int i = 0; i < regions.size() - 1; i++) {
                if (regions.get(i + 1).getFirst() - regions.get(i).getSecond() > novelDistanceSplit) {
                    sections.add(Pair.create(subcontigStart, regions.get(i).getSecond() + window));

                    subcontigStart = regions.get(i + 1).getFirst() - window;
                }
            }

            sections.add(Pair.create(subcontigStart, subcontigStop));

            List<Triple<Integer, Integer, List<CortexVertex>>> wss = new ArrayList<>();
            for (Pair<Integer, Integer> section : sections) {
                List<CortexVertex> ws = new ArrayList<>();
                for (int i = section.getFirst(); i <= section.getSecond(); i++) {
                    ws.add(w.get(i));
                }

                wss.add(Triple.of(section.getFirst(), section.getSecond(), ws));
            }

            return wss;
        }

        return null;
    }

    @NotNull
    private List<Pair<Integer, Integer>> getRegions(Set<CanonicalKmer> rois, List<CortexVertex> cvs) {
        List<Pair<Integer, Integer>> regions = new ArrayList<>();

        int regionStart = -1;
        int regionStop = 0;
        for (int i = 0; i < cvs.size(); i++) {
            CanonicalKmer currentKmer = cvs.get(i).getCanonicalKmer();

            if (rois.contains(currentKmer)) {
                if (regionStart == -1) { regionStart = i; }
                regionStop = i;
            } else {
                if (regionStart > -1) {
                    regions.add(new Pair<>(regionStart, regionStop));

                    regionStart = -1;
                    regionStop = 0;
                }
            }
        }

        if (regionStart > -1) {
            regions.add(new Pair<>(regionStart, regionStop));
        }

        return regions;
    }
}
