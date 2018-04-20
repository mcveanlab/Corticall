package uk.ac.ox.well.cortexjdk.commands.discover.call;

import com.google.common.base.Joiner;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
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

    @Argument(fullName="contigName", shortName="cn", doc="Contigs to process", required=false)
    public HashSet<String> CONTIG_NAMES;

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
                    log.info("  section {}/{}", sectionIndex + 1, sections.size());

                    Triple<Integer, Integer, List<CortexVertex>> section = sections.get(sectionIndex);

                    List<CortexVertex> ws = section.getRight();

                    Map<String, Map<String, String>> allTargets = new HashMap<>();
                    allTargets.put("all", new HashMap<>());
                    for (Set<String> parentName : Arrays.asList(MOTHER, FATHER)) {
                        Map<String, String> parentalTargets = fasterAssembleCandidateHaplotypes(ws, parentName);

                        allTargets.get("all").putAll(parentalTargets);
                        allTargets.put(parentName.iterator().next(), parentalTargets);

                        for (String cn : parentalTargets.keySet()) {
                            log.info("    target background={} name={} length={}", parentName.iterator().next(), cn, parentalTargets.get(cn).length());
                        }
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
                        log.info("\n{}\n{}", makeNoveltyTrack(rois, trimmedQuery.getRight(), lps, true), ma);

                        List<Pair<Integer, Integer>> nrs = getNoveltyRegions(rois, trimmedQuery.getRight(), lps, true);

                        List<VariantContextBuilder> calls = new ArrayList<>();
                        calls.addAll(callSmallBubbles(lps, nrs, rseq.getName().split(" ")[0], section.getLeft() + trimmedQuery.getLeft(), section.getMiddle() + trimmedQuery.getLeft()));
                        calls.addAll(callLargeBubbles(lps, nrs, targets, rseq.getName().split(" ")[0], section.getLeft() + trimmedQuery.getLeft(), section.getMiddle() + trimmedQuery.getLeft()));
                        calls.addAll(callBreakpoints(lps, nrs, rseq.getName().split(" ")[0], section.getLeft() + trimmedQuery.getLeft(), section.getMiddle() + trimmedQuery.getLeft()));

                        List<VariantContextBuilder> merged = mergeBubbles(lps, calls);

                        for (VariantContextBuilder vcb : merged) {
                            vcb.attribute("targets", targets);
                            vcb.attribute("lps", lps);
                        }

                        vcs.addAll(merged);
                    }
                }
            }

            vcs = mergeBreakpoints(seq, vcs);

            for (VariantContextBuilder vcb : vcs) {
                VariantContextBuilder vcbn = assignCoordinates(vcb);

                vcbn.rmAttribute("targets");
                vcbn.rmAttribute("lps");

                VariantContext vc = vcbn.make();

                svcs.add(vc);

                log.info("{}", vc);
            }
            log.info("");
        }

        writeVariants(rois, svcs, vcw);
    }

    private VariantContextBuilder assignCoordinates(VariantContextBuilder vcb) {
        VariantContextBuilder vcbn = new VariantContextBuilder(vcb);

        List<Triple<String, String, Pair<Integer, Integer>>> lps = (ArrayList<Triple<String, String, Pair<Integer, Integer>>>) vcbn.make().getAttribute("lps");

        int start = vcbn.make().getAttributeAsInt("start", 0);
        StringBuilder prevFlankBuilder = new StringBuilder();
        for (int q = start; q >= 0 && getParentalRow(lps, start) == getParentalRow(lps, q); q--) {
            prevFlankBuilder.insert(0, getParentalColumn(lps, q));
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

        int stop = vcbn.make().getAttributeAsInt("stop", 0);
        StringBuilder nextFlankBuilder = new StringBuilder();
        for (int q = stop; q < lps.get(0).getMiddle().length() && getParentalRow(lps, stop) == getParentalRow(lps, q); q++) {
            nextFlankBuilder.append(getParentalColumn(lps, q));
        }

        String nextBack = lps.get(getParentalRow(lps, stop)).getLeft().split(":")[0];
        String nextFlankFw = nextFlankBuilder.toString();
        List<SAMRecord> nextSrs = sortAlignments(nextBack, nextFlankFw);
        SAMRecord nextSr = nextSrs.size() == 0 ? null : nextSrs.get(0);
        if (nextSr != null) {
            vcbn.attribute("nextChrom", nextSr.getContig());
            vcbn.attribute("nextStart", nextSr.getReferencePositionAtReadPosition(1) + 1);
            vcbn.attribute("nextStop", nextSr.getReferencePositionAtReadPosition(nextSr.getReadLength()) + 1);
            vcbn.attribute("nextStrand", nextSr.getReadNegativeStrandFlag() ? "-" : "+");
        }

        SAMRecord sr = null;
        List<SAMRecord> srs = null;
        if (prevSr != null && nextSr != null) {
            sr = prevSr.getStart() < nextSr.getStart() ? prevSr : nextSr;
            srs = prevSr.getStart() < nextSr.getStart() ? prevSrs : nextSrs;
        } else if (prevSr != null) {
            sr = prevSr;
            srs = prevSrs;
        } else if (nextSr != null) {
            sr = nextSr;
            srs = nextSrs;
        }

        if (sr != null) {
            boolean flip = sr.getReadNegativeStrandFlag();
            List<Allele> alleles = vcbn.getAlleles();

            int alignStart = flip ? sr.getEnd() + 1 : sr.getStart() + 1;

            vcbn.chr(sr.getContig());
            vcbn.start(alignStart);
            vcbn.stop(alignStart + vcb.make().getEnd() - vcb.make().getStart());

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

            vcbn.attribute("alt_loci", Joiner.on(";").join(altLoci));
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

    private Set<VariantContextBuilder> mergeBreakpoints(String seq, Set<VariantContextBuilder> callset) {
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

        if (bnds.size() >= 4 && bnds.size() % 2 == 0) {
            for (int i = 0; i <= bnds.size() - 2; i += 2) {
                VariantContext outer0 = bnds.get(i).make();
                VariantContext inner0 = bnds.get(i+1).make();

                List<Triple<String, String, Pair<Integer, Integer>>> lps0 = (ArrayList<Triple<String, String, Pair<Integer, Integer>>>) outer0.getAttribute("lps");
                int q0 = 0;
                StringBuilder kmer0 = new StringBuilder();
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

                for (int j = i + 2; j <= bnds.size() - 2; j += 2) {
                    VariantContext inner1 = bnds.get(j).make();
                    VariantContext outer1 = bnds.get(j+1).make();

                    List<Triple<String, String, Pair<Integer, Integer>>> lps1 = (ArrayList<Triple<String, String, Pair<Integer, Integer>>>) outer1.getAttribute("lps");
                    int q1 = 0;
                    StringBuilder kmer1 = new StringBuilder();
                    for (int q = 1; q < lps1.size(); q++) {
                        if (lps1.get(q).getLeft().equals(outer1.getAttributeAsString("targetName", ""))) {
                            int pos = outer1.getAttributeAsInt("start", 0);
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

                                    String childContig = seq.substring(seq.indexOf(kmer0.toString()) + GRAPH.getKmerSize(), seq.indexOf(kmer1.toString()));

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
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        Set<VariantContextBuilder> newcalls = new LinkedHashSet<>();
        for (VariantContextBuilder vcb : calls) {
            if (!replacements.containsKey(vcb.make().getID())) {
                newcalls.add(vcb);
            } else {
                if (replacements.get(vcb.make().getID()) != null) {
                    newcalls.add(replacements.get(vcb.make().getID()));
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
                                    .attribute("contigName", contigName)
                                    .attribute("prevBase", prevBase)
                                    .attribute("nextBase", nextBase);

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
                                .attribute("contigName", contigName)
                                .attribute("prevBase", prevBase)
                                .attribute("nextBase", nextBase)
                                .attribute("targetName", lps.get(partners.getFirst()).getLeft())
                                .attribute("targetStart", lps.get(partners.getFirst()).getRight().getFirst())
                                .attribute("targetStop", lps.get(partners.getFirst()).getRight().getSecond())
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
                                .attribute("contigName", contigName)
                                .attribute("prevBase", prevBase)
                                .attribute("nextBase", nextBase)
                                .attribute("targetName", lps.get(partners.getSecond()).getLeft())
                                .attribute("targetStart", lps.get(partners.getSecond()).getRight().getFirst())
                                .attribute("targetStop", lps.get(partners.getSecond()).getRight().getSecond())
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

                        List<Allele> alleleStrings = Arrays.asList(Allele.create(pBuilder.toString(), true), Allele.create(cBuilder.toString()));

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
                                .attribute("contigName", contigName)
                                .attribute("prevBase", prevBase)
                                .attribute("nextBase", nextBase);

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
                .rmAttribute("NOVELS")
                .id(id)
                .make()
            );

            for (String sk : vc.getAttributeAsString("NOVELS", "").split(",")) {
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
            if (CONTIG_NAMES == null || CONTIG_NAMES.contains(name[0])) {
                rseqs.add(rseq);
            }
        }
        return rseqs;
    }

    private List<SAMRecord> sortAlignments(String background, String target) {
        List<SAMRecord> a = BACKGROUNDS.get(background).align(target);
        a.removeIf(s -> s.getMappingQuality() == 0);

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

    private String makeNoveltyTrack(Set<CanonicalKmer> rois, String query, List<Triple<String, String, Pair<Integer, Integer>>> lps, boolean expand) {
        int maxLength = 0;
        for (Triple<String, String, Pair<Integer, Integer>> lp : lps) {
            String name = String.format("%s (%d-%d)", lp.getLeft(), lp.getRight().getFirst(), lp.getRight().getSecond());
            maxLength = Math.max(maxLength, name.length());
        }

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

    private List<Pair<Integer, Integer>> getNoveltyRegions(Set<CanonicalKmer> rois, String query, List<Triple<String, String, Pair<Integer, Integer>>> lps, boolean expand) {
        String noveltyTrack = makeNoveltyTrack(rois, query, lps, expand);
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

                if (gs != null) {
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
                if (actualEnd == -1) { actualEnd = w.size() - 1; }

                if (shared > 0) {
                    String contig = TraversalUtils.toContig(w.subList(actualStart, actualEnd));

                    contigs.add(contig);
                }
            }

            int i = 0;
            for (String contig : contigs) {
                String id = String.format("%s:%s_unknown:%s_contig%d_%s", parentName.iterator().next(), parentName.iterator().next(), parentName.iterator().next(), i, "fastasm");

                targets.put(id, contig);

                i++;
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

    /*
    @NotNull
    private Map<String, String> assembleCandidateHaplotypes(Set<CanonicalKmer> rois, List<CortexVertex> ws, List<Pair<Integer, Integer>> regions, Set<String> parentName) {
        Map<String, String> targets = new HashMap<>();

        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = new DirectedWeightedPseudograph<>(CortexEdge.class);

        List<List<CortexVertex>> walks = new ArrayList<>();
        List<String> contigs = new ArrayList<>();

        // First build the non-novel stretches
        for (int i = 0; i < regions.size(); i++) {
            int lastEnd = i == 0 ? 0 : regions.get(i-1).getSecond();
            int nextStart = i == regions.size() - 1 ? ws.size() - 1 : regions.get(i+1).getFirst();

            DirectedWeightedPseudograph<CortexVertex, CortexEdge> gl = assembleLeft(rois, parentName, ws, regions.get(i), lastEnd);
            DirectedWeightedPseudograph<CortexVertex, CortexEdge> gr = assembleRight(rois, parentName, ws, regions.get(i), nextStart);

            Graphs.addGraph(g, gl);
            Graphs.addGraph(g, gr);

            List<CortexVertex> wl = gl.vertexSet().size() == 0 ? new ArrayList<>() : TraversalUtils.toWalk(gl, gl.vertexSet().iterator().next().getKmerAsString(), gl.edgeSet().iterator().next().getColor());
            List<CortexVertex> wr = gr.vertexSet().size() == 0 ? new ArrayList<>() : TraversalUtils.toWalk(gr, gr.vertexSet().iterator().next().getKmerAsString(), gr.edgeSet().iterator().next().getColor());

            if (wl.size() > 0 && !contigs.contains(TraversalUtils.toContig(wl))) {
                walks.add(wl);
                contigs.add(TraversalUtils.toContig(wl));
            }

            if (wr.size() > 0 && !contigs.contains(TraversalUtils.toContig(wr))) {
                walks.add(wr);
                contigs.add(TraversalUtils.toContig(wr));
            }
        }

        targets.putAll(assembleGapFilledCandidates(parentName, g, walks, ws));

        //if (targets.size() == 0) {
        //    targets.putAll(assembleDovetailCandidates(parentName, g, walks));
        //}

        if (targets.size() == 0) {
            for (int i = 0; i < walks.size(); i++) {
                List<CortexVertex> w = walks.get(i);

                String cc = TraversalUtils.toContig(w);
                String id = String.format("%s:%s_unknown:contig%d_%s", parentName.iterator().next(), parentName.iterator().next(), i, "gapped");

                targets.put(id, cc);
            }
        }

        return targets;
    }
    */

    /*
    private Map<String, String> assembleDovetailCandidates(Set<String> parentName, DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, List<List<CortexVertex>> walks) {
        // Now try dovetail overlaps with remaining gaps
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> gm = new DirectedWeightedPseudograph<>(CortexEdge.class);
        Graphs.addGraph(gm, g);

        int numDovetailOverlaps = 0;

        for (int i = 0; i < walks.size() - 1; i++) {
            String cb = TraversalUtils.toContig(walks.get(i));
            String ca = TraversalUtils.toContig(walks.get(i+1));

            int longestMatch = 0;
            for (int j = 1; j < cb.length() && j < ca.length(); j++) {
                if (cb.substring(cb.length() - j, cb.length()).equals(ca.substring(0, j))) {
                    longestMatch = j;
                }
            }

            if (longestMatch > 11) {
                String cc = cb + ca.substring(longestMatch, ca.length());

                CortexVertex v0 = null;
                for (int j = 0; j <= cc.length() - GRAPH.getKmerSize(); j++) {
                    String sk = cc.substring(j, j + GRAPH.getKmerSize());
                    CortexRecord cr = GRAPH.findRecord(sk);

                    CortexVertex v = new CortexVertexFactory()
                            .bases(sk)
                            .record(cr)
                            .make();

                    gm.addVertex(v);

                    if (v0 != null) {
                        gm.addEdge(v0, v, new CortexEdge(v0, v, GRAPH.getColorsForSampleNames(parentName).iterator().next(), 1.0f));
                    }

                    v0 = v;
                }

                numDovetailOverlaps++;
            }
        }

        Map<String, String> targets = new HashMap<>();
        if (numDovetailOverlaps > 0) {
            targets.putAll(extractTargets(parentName, gm, "dovetail"));
        }

        return targets;
    }
    */

    /*
    private Map<String, String> assembleGapFilledCandidates(Set<String> parentName, DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, List<List<CortexVertex>> walks, List<CortexVertex> ws) {
        DirectedGraph<CortexByteKmer, DefaultEdge> gw = new DefaultDirectedGraph<>(DefaultEdge.class);

        gw.addVertex(ws.get(0).getKmerAsByteKmer());
        for (int i = 1; i < ws.size(); i++) {
            gw.addVertex(ws.get(i).getKmerAsByteKmer());
            gw.addEdge(ws.get(i-1).getKmerAsByteKmer(), ws.get(i).getKmerAsByteKmer());
        }

        // Now bridge any gaps that we can
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> gm = new DirectedWeightedPseudograph<>(CortexEdge.class);
        Graphs.addGraph(gm, g);

        int numGapsFilled = 0;

        for (int i = 0; i < walks.size() - 1; i++) {
            DirectedWeightedPseudograph<CortexVertex, CortexEdge> gg = new DirectedWeightedPseudograph<>(CortexEdge.class);
            for (int j = i + 1; j < walks.size() && gg.vertexSet().size() == 0; j++) {
                gg = assembleMiddle(parentName, walks.get(i), walks.get(j), gw);
            }

            if (gg.vertexSet().size() > 0) {
                Graphs.addGraph(gm, gg);
                numGapsFilled++;
            } else {
                TraversalEngine ef = new TraversalEngineFactory()
                        .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                        .traversalDirection(FORWARD)
                        .combinationOperator(OR)
                        .stoppingRule(ContigStopper.class)
                        .maxBranchLength(500)
                        .graph(GRAPH)
                        .links(LINKS)
                        .make();

                TraversalEngine er = new TraversalEngineFactory()
                        .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                        .traversalDirection(REVERSE)
                        .combinationOperator(OR)
                        .stoppingRule(ContigStopper.class)
                        .maxBranchLength(500)
                        .graph(GRAPH)
                        .links(LINKS)
                        .make();

                if (i < walks.size() - 1) {
                    Graphs.addGraph(gm, ef.dfs(walks.get(i).get(walks.get(i).size() - 1).getKmerAsString()));
                }

                if (i > 0) {
                    Graphs.addGraph(gm, er.dfs(walks.get(i).get(walks.get(i).size() - 1).getKmerAsString()));
                }
            }
        }

        Map<String, String> targets = new HashMap<>();
        if (numGapsFilled > 0) {
            targets.putAll(extractTargets(parentName, gm, "gapfilled"));
        }

        return targets;
    }
    */

    /*
    private Map<String, String> extractTargets(Set<String> parentName, DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, String suffix) {
        Map<String, String> targets = new HashMap<>();
        ConnectivityInspector<CortexVertex, CortexEdge> ci = new ConnectivityInspector<>(g);
        int targetIndex = 0;
        for (Set<CortexVertex> cs : ci.connectedSets()) {
            List<CortexVertex> cl = TraversalUtils.toWalk(g, cs.iterator().next().getKmerAsString(), g.edgeSet().iterator().next().getColor());

            String cc = TraversalUtils.toContig(cl);
            String id = String.format("%s:%s_unknown:contig%d_%s", parentName.iterator().next(), parentName.iterator().next(), targetIndex, suffix);

            if (cc.length() > 0) {
                targets.put(id, cc);
                targetIndex++;
            }
        }

        return targets;
    }
    */

    /*
    private DirectedWeightedPseudograph<CortexVertex, CortexEdge> assembleMiddle(Set<String> parentName, List<CortexVertex> wLeft, List<CortexVertex> wRight, DirectedGraph<CortexByteKmer, DefaultEdge> gw) {
        int childColor = GRAPH.getColorForSampleName(ROIS.getSampleName(0));

        for (int i = wLeft.size() - 1; i > 0; i--) {
            Map<Integer, Set<CortexByteKmer>> cbkOut = TraversalUtils.getAllNextKmers(wLeft.get(i).getCortexRecord(), !wLeft.get(i).getKmerAsString().equals(wLeft.get(i).getCanonicalKmer().getKmerAsString()));
            Set<CortexByteKmer> outgoingEdges = new HashSet<>();
            for (int c : GRAPH.getColorsForSampleNames(parentName)) {
                outgoingEdges.addAll(cbkOut.get(c));
            }

            if (gw.containsVertex(wLeft.get(i).getKmerAsByteKmer())) {
                outgoingEdges.removeAll(Graphs.successorListOf(gw, wLeft.get(i).getKmerAsByteKmer()));
            } else {
                outgoingEdges.removeAll(cbkOut.get(childColor));
            }

            if (outgoingEdges.size() > 0) {
                for (int j = 0; j < wRight.size() - 1; j++) {
                    Map<Integer, Set<CortexByteKmer>> cbkIn = TraversalUtils.getAllPrevKmers(wRight.get(j).getCortexRecord(), !wRight.get(j).getKmerAsString().equals(wRight.get(j).getCanonicalKmer().getKmerAsString()));
                    Set<CortexByteKmer> incomingEdges = new HashSet<>();
                    for (int c : GRAPH.getColorsForSampleNames(parentName)) {
                        incomingEdges.addAll(cbkIn.get(c));
                    }

                    if (gw.containsVertex(wRight.get(j).getKmerAsByteKmer())) {
                        incomingEdges.removeAll(Graphs.predecessorListOf(gw, wRight.get(j).getKmerAsByteKmer()));
                    } else {
                        incomingEdges.removeAll(cbkIn.get(childColor));
                    }

                    if (incomingEdges.size() > 0) {
                        TraversalEngine ef = new TraversalEngineFactory()
                                .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                                .traversalDirection(FORWARD)
                                .combinationOperator(OR)
                                .stoppingRule(DestinationStopper.class)
                                .maxBranchLength(2000)
                                .graph(GRAPH)
                                .links(LINKS)
                                .make();

                        TraversalEngine er = new TraversalEngineFactory()
                                .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                                .traversalDirection(REVERSE)
                                .combinationOperator(OR)
                                .stoppingRule(DestinationStopper.class)
                                .maxBranchLength(2000)
                                .graph(GRAPH)
                                .links(LINKS)
                                .make();

                        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = ef.dfs(wLeft.get(i-1).getKmerAsString(), wRight.get(j+1).getKmerAsString(), SequenceUtils.reverseComplement(wRight.get(wRight.size() - j - 1).getKmerAsString()));
                        if (g == null || g.vertexSet().size() == 0) {
                            g = er.dfs(wRight.get(j+1).getKmerAsString(), wLeft.get(i-1).getKmerAsString(), SequenceUtils.reverseComplement(wLeft.get(wLeft.size() - i).getKmerAsString()));
                        }

                        if (g != null && g.vertexSet().size() > 0) {
                            return g;
                        }
                    }
                }
            }
        }

        return new DirectedWeightedPseudograph<>(CortexEdge.class);
    }
    */

    /*
    private DirectedWeightedPseudograph<CortexVertex, CortexEdge> assembleLeft(Set<CanonicalKmer> rois, Set<String> parentName, List<CortexVertex> ws, Pair<Integer, Integer> region, int lastEnd) {
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> gl = new DirectedWeightedPseudograph<>(CortexEdge.class);

        for (int j = region.getFirst() - 1; j > lastEnd && !rois.contains(ws.get(j).getCanonicalKmer()); j--) {
            boolean hasCoverage = false;
            for (int c : GRAPH.getColorsForSampleNames(parentName)) {
                if (ws.get(j).getCortexRecord().getCoverage(c) > 0) {
                    hasCoverage = true;
                    break;
                }
            }

            if (hasCoverage) {
                TraversalEngine e = new TraversalEngineFactory()
                        .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                        .traversalDirection(REVERSE)
                        .combinationOperator(OR)
                        .stoppingRule(ContigStopper.class)
                        .maxBranchLength((region.getFirst() - 1) - (lastEnd + 1))
                        .graph(GRAPH)
                        .links(LINKS)
                        .make();

                DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = e.dfs(ws.get(j).getKmerAsString());

                if (g != null) {
                    if (g.vertexSet().size() > gl.vertexSet().size()) {
                        gl = g;
                    } else if (gl.vertexSet().size() > 0 && g.vertexSet().size() <= gl.vertexSet().size()) {
                        break;
                    }
                }
            }
        }

        return gl;
    }
    */

    /*
    private DirectedWeightedPseudograph<CortexVertex, CortexEdge> assembleRight(Set<CanonicalKmer> rois, Set<String> parentName, List<CortexVertex> ws, Pair<Integer, Integer> region, int nextStart) {
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> gr = new DirectedWeightedPseudograph<>(CortexEdge.class);

        for (int j = region.getSecond() + 1; j < nextStart && !rois.contains(ws.get(j).getCanonicalKmer()); j++) {
            boolean hasCoverage = false;
            for (int c : GRAPH.getColorsForSampleNames(parentName)) {
                if (ws.get(j).getCortexRecord().getCoverage(c) > 0) {
                    hasCoverage = true;
                    break;
                }
            }

            if (hasCoverage) {
                TraversalEngine e = new TraversalEngineFactory()
                        .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                        .traversalDirection(FORWARD)
                        .combinationOperator(OR)
                        .stoppingRule(ContigStopper.class)
                        .maxBranchLength((nextStart - 1) - (region.getSecond() + 1))
                        .graph(GRAPH)
                        .links(LINKS)
                        .make();

                DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = e.dfs(ws.get(j).getKmerAsString());

                if (g != null) {
                    if (g.vertexSet().size() > gr.vertexSet().size()) {
                        gr = g;
                    } else if (gr.vertexSet().size() > 0 && g.vertexSet().size() <= gr.vertexSet().size()) {
                        break;
                    }
                }
            }
        }

        return gr;
    }
    */

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
