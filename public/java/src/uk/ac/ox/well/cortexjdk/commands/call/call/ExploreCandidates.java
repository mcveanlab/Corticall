package uk.ac.ox.well.cortexjdk.commands.call.call;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.alignment.sw.SmithWaterman;
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

    @Argument(fullName = "background", shortName = "b", doc = "Background name")
    public HashSet<String> BACKGROUNDS;

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

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing novel kmers...")
                .message("records processed")
                .maxRecord(used.size())
                .make(log);

        int numContigs = 0;
        for (CanonicalKmer ck : used.keySet()) {
            pm.update();

            if (used.get(ck) == null) {
                DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = e.dfs(ck);
                List<CortexVertex> w = TraversalEngine.toWalk(g, ck);

                Pair<Integer, Integer> numMarked = markUsedRois(used, w);

                log.info("    {} seed={} dfs={} contig={} numNewlyMarked={} numAlreadyMarked={}", numContigs, ck, g.vertexSet().size(), w.size(), numMarked.getFirst(), numMarked.getSecond());

                numContigs++;

                int q0 = -1, q1 = -1;
                for (int i = 0; i < w.size(); i++) {
                    if (used.containsKey(w.get(i).getCanonicalKmer())) {
                        if (q0 == -1) { q0 = i; }
                        q1 = i;
                    }
                }

                int s0 = q0 - 100 >= 0 ? q0 - 100 : 0;
                int s1 = q1 + 100 < w.size() - 1 ? q1 + 100 : w.size() - 1;

                List<Integer> colors = new ArrayList<>();
                colors.add(getTraversalColor(GRAPH, ROIS));
                colors.add(getTraversalColor(GRAPH, ROIS));
                colors.addAll(getParentalColors(GRAPH, BACKGROUNDS));

                boolean processedNovelTrack = false;

                for (int c : colors) {
                    Map<Integer, Set<String>> divOutVertices = new TreeMap<>();
                    Map<Integer, Set<String>> divInVertices = new TreeMap<>();

                    StringBuilder sb = new StringBuilder();
                    for (int i = 0; i < w.size(); i++) {
                        CortexVertex v = w.get(i);

                        Set<String> divIn = divergenceDegree(v, false, getTraversalColor(GRAPH, ROIS), c);
                        Set<String> divOut = divergenceDegree(v, true, getTraversalColor(GRAPH, ROIS), c);

                        if (!processedNovelTrack && c == getTraversalColor(GRAPH, ROIS) && used.containsKey(v.getCanonicalKmer())) {
                            sb.append("!");
                        } else if (v.getCortexRecord().getCoverage(c) == 0) {
                            sb.append(" ");
                        } else if (divIn.size() + divOut.size() > 0) {
                            if (divIn.size() > 0 && divOut.size() == 0) {
                                sb.append("\\");
                                if (i >= s0 && i <= s1) {
                                    divInVertices.put(i, divIn);
                                }
                            } else if (divIn.size() == 0 && divOut.size() > 0) {
                                sb.append("/");
                                if (i >= s0 && i <= s1) {
                                    divOutVertices.put(i, divOut);
                                }
                            } else {
                                sb.append("=");
                                if (i >= s0 && i <= s1) {
                                    divInVertices.put(i, divIn);
                                    divOutVertices.put(i, divOut);
                                }
                            }
                        } else {
                            sb.append("_");
                        }
                    }

                    if (c == getTraversalColor(GRAPH, ROIS)) { processedNovelTrack = true; }

                    String s = sb.toString();

                    log.info("{} {}: {}{}{}", String.format("%-3d", c), String.format("%-20s", GRAPH.getSampleName(c)), s0 == 0 ? "" : "...", s.substring(s0, s1), s1 == s.length() ? "" : "...");

                    if (divOutVertices.size() > 0 && divInVertices.size() > 0) {
                        for (int ovi : divOutVertices.keySet()) {
                            Set<String> ovs = divOutVertices.get(ovi);
                            Set<String> ivs = new HashSet<>();

                            for (int ivi : divInVertices.keySet()) {
                                if (ivi > ovi) {
                                    ivs.addAll(divInVertices.get(ivi));
                                }
                            }

                            TraversalEngine eo = new TraversalEngineFactory()
                                    .traversalColor(c)
                                    .traversalDirection(FORWARD)
                                    .combinationOperator(OR)
                                    .graph(GRAPH)
                                    .links(LINKS)
                                    //.references(REFERENCES)
                                    .rois(ROIS)
                                    .sink(ivs)
                                    .stoppingRule(DestinationStopper.class)
                                    .make();

                            for (String ov : ovs) {
                                DirectedWeightedPseudograph<CortexVertex, CortexEdge> go = eo.dfs(ov);
                                List<CortexVertex> gw = TraversalEngine.toWalk(go, ov);
                                String gc = TraversalEngine.toContig(gw);

                                if (go != null) {
                                    for (CortexVertex cv : go.vertexSet()) {
                                        if (ivs.contains(cv.getKmerAsString())) {
                                            for (int ivi : divInVertices.keySet()) {
                                                List<CortexVertex> subContig = new ArrayList<>();
                                                for (int i = ovi + 1; i < ivi; i++) {
                                                    subContig.add(w.get(i));
                                                }

                                                String sc = TraversalEngine.toContig(subContig);

                                                if (divInVertices.get(ivi).contains(cv.getKmerAsString())) {
                                                    log.info("{}{} {} {}", String.format("%-26s", ""), ovi, ivi, sc);
                                                    log.info("{}{} {} {}", String.format("%-26s", ""), ovi, ivi, gc);

                                                    trimToAlleles(gc, sc);
                                                }
                                            }
                                        }
                                    }

                                    /*
                                    List<CortexVertex> subContig = new ArrayList<>();
                                    for (int i = ovi + 1; i < ivi; i++) {
                                        subContig.add(w.get(i));
                                    }
                                    log.info("alt {} {} {}{}", ovi, ivi, String.format("%-26s", ""), TraversalEngine.toContig(subContig));
                                    log.info("hap {} {} {}{}", ovi, ivi, String.format("%-26s", ""), TraversalEngine.toContig(TraversalEngine.toWalk(go, ov)));
                                    */
                                }
                            }
                        }
                    }
                }

                log.info("");
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

    private void trimToAlleles(String altContig, String compContig) {
        SmithWaterman sw = new SmithWaterman(compContig, altContig);
        String[] aligns = sw.getAlignment();

        List<VariantContext> vcs = new ArrayList<>();

        int alleleStart = 0, alleleStop = 0;
        StringBuilder altBuilder = new StringBuilder();
        StringBuilder compBuilder = new StringBuilder();
        for (int i = 0; i < aligns[0].length(); i++, alleleStop++) {
            if (aligns[0].charAt(i) != aligns[1].charAt(i)) {
                if (altBuilder.length() == 0) {
                    alleleStart = i;
                }

                altBuilder.append(aligns[0].charAt(i));
                compBuilder.append(aligns[1].charAt(i));
            } else {
                // add variant to pile
                if (altBuilder.length() > 0) {
                    log.info("{} {} {}", alleleStart, altBuilder, compBuilder);

                    Set<Character> a = new HashSet<>();
                    for (int q = 0; q < altBuilder.length(); q++) { a.add(altBuilder.charAt(q)); }

                    Set<Character> b = new HashSet<>();
                    for (int q = 0; q < compBuilder.length(); q++) { b.add(compBuilder.charAt(q)); }

                    if (a.contains('-') && (a.contains('A') || a.contains('C') || a.contains('G') || a.contains('T'))) {
                        log.info("{}", a);
                    }

                    if (b.contains('-') && (b.contains('A') || b.contains('C') || b.contains('G') || b.contains('T'))) {
                        log.info("{}", b);
                    }

                    vcs.add(buildVariantContext(aligns, alleleStart, altBuilder, compBuilder));

                    altBuilder = new StringBuilder();
                    compBuilder = new StringBuilder();
                }
            }
        }

        if (altBuilder.length() > 0) {
            vcs.add(buildVariantContext(aligns, alleleStart, altBuilder, compBuilder));
        }
    }

    private VariantContext buildVariantContext(String[] aligns, int alleleStart, StringBuilder altBuilder, StringBuilder compBuilder) {
        String altAllele = altBuilder.toString();
        String compAllele = compBuilder.toString();

        if (altAllele.contains("-")) {
            altBuilder = new StringBuilder(altAllele.replaceAll("-", ""));
            altBuilder.insert(0, aligns[0].charAt(alleleStart - 1));
            altAllele = altBuilder.toString();

            compBuilder.insert(0, aligns[1].charAt(alleleStart - 1));
            compAllele = compBuilder.toString();
        } else if (compAllele.contains("-")) {
            compBuilder = new StringBuilder(compAllele.replaceAll("-", ""));
            compBuilder.insert(0, aligns[1].charAt(alleleStart - 1));
            compAllele = compBuilder.toString();

            altBuilder.insert(0, aligns[0].charAt(alleleStart - 1));
            altAllele = altBuilder.toString();
        }

        VariantContext vc = new VariantContextBuilder()
                .chr("unknown")
                .start(alleleStart)
                .computeEndFromAlleles(Arrays.asList(Allele.create(compAllele, true), Allele.create(altAllele, false)), alleleStart)
                .alleles(compAllele, altAllele)
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

    private Collection<Integer> getParentalColors(DeBruijnGraph graph, Set<String> parentNames) {
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
