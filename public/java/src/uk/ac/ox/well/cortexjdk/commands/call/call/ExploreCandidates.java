package uk.ac.ox.well.cortexjdk.commands.call.call;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.sw.SWResult;
import uk.ac.ox.well.cortexjdk.utils.alignment.swold.SmithWaterman;
import uk.ac.ox.well.cortexjdk.utils.alignment.sw.SmithWaterman2;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.DestinationStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.NovelContinuationStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.*;

import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;

public class ExploreCandidates extends Module {
    @Argument(fullName = "graph", shortName = "g", doc = "Graph")
    public CortexGraph GRAPH;

    @Argument(fullName = "links", shortName = "l", doc = "Links", required=false)
    public ArrayList<CortexLinks> LINKS;

    //@Argument(fullName = "references", shortName = "R", doc = "References", required=false)
    //public ArrayList<IndexedReference> REFERENCES;

    @Argument(fullName = "roi", shortName = "r", doc = "ROI")
    public CortexGraph ROIS;

    @Argument(fullName = "background", shortName = "b", doc = "Background color name")
    public HashSet<String> BACKGROUNDS;

    @Argument(fullName = "windowAroundRois", shortName = "w", doc = "Window of consideration around variants")
    public Integer WINDOW = 100;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(getTraversalColor(GRAPH, ROIS))
                .traversalDirection(BOTH)
                .combinationOperator(OR)
                .graph(GRAPH)
                .links(LINKS)
                //.references(REFERENCES)
                .rois(ROIS)
                .stoppingRule(NovelContinuationStopper.class)
                .make();

        Map<CanonicalKmer, List<CortexVertex>> used = loadRois(ROIS);

        List<Integer> colors = new ArrayList<>();
        colors.add(getTraversalColor(GRAPH, ROIS));
        colors.addAll(getBackgroundColors(GRAPH, BACKGROUNDS));

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing novel kmers...")
                .message("records processed")
                .maxRecord(used.size())
                .make(log);

        int numContigs = 0;
        for (CanonicalKmer ck : used.keySet()) {
            pm.update();

            if (used.get(ck) == null) {
                List<CortexVertex> w = e.walk(ck);

                Pair<Integer, Integer> numMarked = markUsedRois(used, w);
                Pair<Integer, Integer> range = getNovelKmerRange(used, w);

                log.info("  * contig={} seed={} len={} numNewlyMarked={} numAlreadyMarked={}", numContigs, ck, w.size(), numMarked.getFirst(), numMarked.getSecond());

                for (int c : colors) {
                    displayContigAndAnnotations(w, range, c, used);

                    Map<Integer, VariantContext> nvcs = callVariantsAgainstBackground(w, range, c, used);

                    log.info("        backgroundColor={} backgroundName={} nvcs={}", c, GRAPH.getSampleName(c), nvcs.size());

                    for (int offset : nvcs.keySet()) {
                        log.info("        w={} v0={} v1={} offset={} vc={}", w.size(), range.getFirst(), range.getSecond(), offset, nvcs.get(offset));
                    }

                    log.info("    ==");
                    log.info("");
                }

                numContigs++;
            }
        }

        int numNovelKmersAssigned = 0;
        for (CanonicalKmer ck : used.keySet()) {
            if (used.get(ck) != null) {
                numNovelKmersAssigned++;
            }
        }

        log.info("Assigned {}/{} novel kmers to {} contigs", numNovelKmersAssigned, used.size(), numContigs);
    }

    private void displayContigAndAnnotations(List<CortexVertex> w, Pair<Integer, Integer> range, int c, Map<CanonicalKmer, List<CortexVertex>> used) {
        StringBuilder t0 = new StringBuilder();
        StringBuilder t1 = new StringBuilder();
        StringBuilder t2 = new StringBuilder();
        List<CortexVertex> wsub = new ArrayList<>();
        for (int i = range.getFirst(); i <= range.getSecond(); i++) {
            wsub.add(w.get(i));

            t0.append(w.get(i).getKmerAsString().substring(0, 1));
            t1.append(used.containsKey(w.get(i).getCanonicalKmer()) ? "!" : "_");
            int divOut = divergenceDegree(w.get(i), true, getTraversalColor(GRAPH, ROIS), c).size();
            int divIn = divergenceDegree(w.get(i), false, getTraversalColor(GRAPH, ROIS), c).size();

            if      (divOut >  0 && divIn == 0) { t2.append("\\");  }
            else if (divOut == 0 && divIn >  0) { t2.append("/"); }
            else if (divOut >  0 && divIn >  0) { t2.append("=");  }
            else                                { t2.append(" ");  }
        }

        log.info("    {}", t0.toString());
        log.info("    {}", t1.toString());
        log.info("    {}", t2.toString());
    }

    private Map<Integer, VariantContext> callVariantsAgainstBackground(List<CortexVertex> w, Pair<Integer, Integer> range, int c, Map<CanonicalKmer, List<CortexVertex>> used) {
        List<VariantContext> vcs = new ArrayList<>();

        Map<Integer, Set<String>> divOutVertices = new TreeMap<>();
        Map<Integer, Set<String>> divInVertices = new TreeMap<>();

        for (int i = 0; i < w.size(); i++) {
            CortexVertex v = w.get(i);

            if (i >= range.getFirst() && i <= range.getSecond()) {
                Set<String> divIn = divergenceDegree(v, false, getTraversalColor(GRAPH, ROIS), c);
                Set<String> divOut = divergenceDegree(v, true, getTraversalColor(GRAPH, ROIS), c);

                if (divIn.size() + divOut.size() > 0) {
                    if (divIn.size() > 0) { divInVertices.put(i, divIn); }
                    if (divOut.size() > 0) { divOutVertices.put(i, divOut); }
                }
            }
        }

        if (divOutVertices.size() > 0 && divInVertices.size() > 0) {
            for (int ovi : divOutVertices.keySet()) {
                Set<String> ovs = divOutVertices.get(ovi);
                Set<String> ivs = new HashSet<>();

                for (int ivi : divInVertices.keySet()) {
                    if (ivi > ovi) {
                        ivs.addAll(divInVertices.get(ivi));
                    }
                }

                for (String ov : ovs) {
                    int numNovelsFinal = 0;
                    List<CortexVertex> traversalWalkFinal = null;
                    List<CortexVertex> backgroundWalkFinal = null;
                    String ivFinal = null;

                    for (String iv : ivs) {
                        TraversalEngine em = new TraversalEngineFactory()
                                .traversalColor(c)
                                .traversalDirection(FORWARD)
                                .combinationOperator(OR)
                                .graph(GRAPH)
                                .links(LINKS)
                                .rois(ROIS)
                                .sink(iv)
                                .stoppingRule(DestinationStopper.class)
                                .make();

                        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = em.dfs(ov);

                        if (g != null) {
                            List<CortexVertex> backgroundWalk = TraversalEngine.toWalk(g, ov);

                            int numNovels = 0;
                            List<CortexVertex> traversalWalk = new ArrayList<>();
                            for (int ivi : divInVertices.keySet()) {
                                if (divInVertices.get(ivi).contains(iv)) {
                                    for (int j = ovi + 1; j < ivi; j++) {
                                        if (used.containsKey(w.get(j).getCanonicalKmer())) {
                                            numNovels++;
                                        }

                                        traversalWalk.add(w.get(j));
                                    }

                                    break;
                                }
                            }

                            if (backgroundWalkFinal == null || (numNovels >= numNovelsFinal && backgroundWalk.size() <= backgroundWalkFinal.size())) {
                                numNovelsFinal = numNovels;
                                traversalWalkFinal = traversalWalk;
                                backgroundWalkFinal = backgroundWalk;
                                ivFinal = iv;
                            }
                        }
                    }

                    if (traversalWalkFinal != null && backgroundWalkFinal != null && numNovelsFinal > 0) {
                        String backgroundColorContig = TraversalEngine.toContig(backgroundWalkFinal);
                        String traversalColorContig = TraversalEngine.toContig(traversalWalkFinal);

                        vcs.addAll(trimToAlleles(traversalColorContig, backgroundColorContig, ovi + GRAPH.getKmerSize() + 1, ov, ivFinal, numNovelsFinal));
                    }
                }
            }
        }

        Map<Integer, VariantContext> nvcs = new TreeMap<>();
        for (VariantContext vc : vcs) {
            nvcs.put(vc.getStart() - range.getFirst() - GRAPH.getKmerSize(), vc);
        }

        return nvcs;
    }

    private Pair<Integer, Integer> getNovelKmerRange(Map<CanonicalKmer, List<CortexVertex>> used, List<CortexVertex> w) {
        int q0 = -1, q1 = -1;
        for (int i = 0; i < w.size(); i++) {
            if (used.containsKey(w.get(i).getCanonicalKmer())) {
                if (q0 == -1) { q0 = i; }
                q1 = i;
            }
        }

        int s0 = q0 - WINDOW >= 0 ? q0 - WINDOW : 0;
        int s1 = q1 + WINDOW < w.size() - 1 ? q1 + WINDOW : w.size() - 1;

        return new Pair<>(s0, s1);
    }

    private List<VariantContext> trimToAlleles(String traversalColorContig, String backgroundColorContig, int offset, String ov, String iv, int numNovels) {
        SmithWaterman2 sw = new SmithWaterman2();
        String[] aligns = sw.getAlignment(traversalColorContig, backgroundColorContig);

        //log.info("{}", traversalColorContig);
        //log.info("{}", backgroundColorContig);
        //log.info("{}", aligns[0]);
        //log.info("{}", aligns[1]);

        List<VariantContext> vcs = new ArrayList<>();

        int currentPos = 0;
        int alleleStart = -1;
        StringBuilder altBuilder = new StringBuilder();
        StringBuilder compBuilder = new StringBuilder();
        for (int i = 0; i < aligns[0].length(); i++) {
            if (aligns[0].charAt(i) != aligns[1].charAt(i)) {
                if (altBuilder.length() == 0) {
                    alleleStart = currentPos;
                }

                altBuilder.append(aligns[0].charAt(i));
                compBuilder.append(aligns[1].charAt(i));
            } else {
                if (altBuilder.length() > 0) {
                    vcs.add(buildVariantContext(aligns, alleleStart, altBuilder, compBuilder, offset, ov, iv, numNovels));

                    alleleStart = -1;
                    altBuilder = new StringBuilder();
                    compBuilder = new StringBuilder();
                }
            }

            if (aligns[0].charAt(i) != '-') {
                currentPos++;
            }
        }

        if (altBuilder.length() > 0) {
            vcs.add(buildVariantContext(aligns, alleleStart, altBuilder, compBuilder, offset, ov, iv, numNovels));
        }

        return vcs;
    }

    private VariantContext buildVariantContext(String[] aligns, int alleleStart, StringBuilder altBuilder, StringBuilder compBuilder, int offset, String ov, String iv, int numNovels) {
        String altAllele = altBuilder.toString().replaceAll("-", "");
        String compAllele = compBuilder.toString().replaceAll("-", "");

        if (altAllele.length() != compAllele.length()) {
            alleleStart -= 1;

            altBuilder = new StringBuilder(altAllele);
            altBuilder.insert(0, aligns[0].replaceAll("-", "").charAt(alleleStart));
            altAllele = altBuilder.toString();

            compBuilder = new StringBuilder(compAllele);
            compBuilder.insert(0, aligns[0].replaceAll("-", "").charAt(alleleStart));
            compAllele = compBuilder.toString();
        }

        VariantContext vc = new VariantContextBuilder()
                .chr("unknown")
                .start(alleleStart + offset)
                .computeEndFromAlleles(Arrays.asList(Allele.create(altAllele, true), Allele.create(compAllele, false)), alleleStart + offset)
                .alleles(altAllele, compAllele)
                .attribute("ov", ov)
                .attribute("iv", iv)
                .attribute("numNovels", numNovels)
                .make();

        return vc;
    }

    private Pair<Integer, Integer> markUsedRois(Map<CanonicalKmer, List<CortexVertex>> used, List<CortexVertex> w) {
        int numNewlyMarked = 0, numAlreadyMarked = 0;
        for (CortexVertex v : w) {
            if (used.containsKey(v.getCanonicalKmer()) && used.get(v.getCanonicalKmer()) != null) {
                numAlreadyMarked++;
            }
        }

        for (CortexVertex v : w) {
            if (used.containsKey(v.getCanonicalKmer())) {
                if ((used.get(v.getCanonicalKmer()) == null || w.size() > used.get(v.getCanonicalKmer()).size())) {
                    used.put(v.getCanonicalKmer(), w);
                    numNewlyMarked++;
                }
            }
        }

        return new Pair<>(numNewlyMarked, numAlreadyMarked);
    }

    private Map<CanonicalKmer, List<CortexVertex>> loadRois(CortexGraph rois) {
        Map<CanonicalKmer, List<CortexVertex>> used = new HashMap<>();
        for (CortexRecord cr : rois) {
            used.put(cr.getCanonicalKmer(), null);
        }

        return used;
    }

    private int getTraversalColor(DeBruijnGraph graph, CortexGraph rois) {
        return graph.getColorForSampleName(rois.getSampleName(0));
    }

    private Collection<Integer> getBackgroundColors(DeBruijnGraph graph, Set<String> parentNames) {
        return graph.getColorsForSampleNames(parentNames);
    }

    private Set<String> divergenceDegree(CortexVertex v, boolean goForward, int traversalColor, int comparisonColor) {
        Map<Integer, Set<CortexByteKmer>> outEdges = TraversalEngine.getAllNextKmers(v.getCortexRecord(), !v.getKmerAsString().equalsIgnoreCase(v.getCanonicalKmer().getKmerAsString()));
        Map<Integer, Set<CortexByteKmer>> inEdges = TraversalEngine.getAllPrevKmers(v.getCortexRecord(), !v.getKmerAsString().equalsIgnoreCase(v.getCanonicalKmer().getKmerAsString()));
        Map<Integer, Set<CortexByteKmer>> edges = goForward ? outEdges : inEdges;

        Set<CortexByteKmer> remainingEdges = edges.get(comparisonColor);
        remainingEdges.removeAll(edges.get(traversalColor));

        Set<String> remainingEdgeStrs = new HashSet<>();
        for (CortexByteKmer cbk : remainingEdges) {
            remainingEdgeStrs.add(new String(cbk.getKmer()));
        }

        return remainingEdgeStrs;
    }
}
