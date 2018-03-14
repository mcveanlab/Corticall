package uk.ac.ox.well.cortexjdk.commands.discover.call;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;
import org.jgrapht.Graphs;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.mapdb.Atomic;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.mosaic.MosaicAligner;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.containers.ContainerUtils;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContigStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.DestinationStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.*;

import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.REVERSE;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalUtils.getAllNextKmers;

public class FindPaths extends Module {
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

    @Argument(fullName="reference", shortName="R", doc="Reference", required=false)
    public HashMap<String, IndexedReference> REFERENCE;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Set<CanonicalKmer> rois = loadRois(ROIS);

        ReferenceSequence rseq;
        while ((rseq = PARTITIONS.nextSequence()) != null) {
            log.info("{}", rseq.getName());

            List<CortexVertex> w = loadChildWalk(rseq, GRAPH, LINKS);
            List<Triple<Integer, Integer, List<CortexVertex>>> sections = sectionContig(rois, w, 100, 500);

            for (Triple<Integer, Integer, List<CortexVertex>> section : sections) {
                List<CortexVertex> ws = section.getRight();
                String query = TraversalUtils.toContig(ws);

                Map<String, String> targets = new HashMap<>();
                List<Pair<Integer, Integer>> regions = getRegions(rois, ws);

                for (Set<String> parentName : Arrays.asList(MOTHER, FATHER)) {
                    DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = new DirectedWeightedPseudograph<>(CortexEdge.class);
                    List<List<CortexVertex>> walks = new ArrayList<>();
                    List<String> contigs = new ArrayList<>();

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

                    for (int i = 0; i < walks.size() - 1; i++) {
                        DirectedWeightedPseudograph<CortexVertex, CortexEdge> gm = new DirectedWeightedPseudograph<>(CortexEdge.class);
                        for (int j = i + 1; j < walks.size() && gm.vertexSet().size() == 0; j++) {
                            gm = assembleMiddle(parentName, walks.get(i), walks.get(j));
                        }

                        Graphs.addGraph(g, gm);

                        TraversalEngine e = new TraversalEngineFactory()
                                .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                                .traversalDirection(BOTH)
                                .combinationOperator(OR)
                                .stoppingRule(ContigStopper.class)
                                .maxBranchLength(2*GRAPH.getKmerSize())
                                .graph(GRAPH)
                                .links(LINKS)
                                .make();

                        DirectedWeightedPseudograph<CortexVertex, CortexEdge> gb = e.dfs(walks.get(i).get(0).getKmerAsString());
                        if (gb != null) { Graphs.addGraph(g, gb); }

                        DirectedWeightedPseudograph<CortexVertex, CortexEdge> gf = e.dfs(walks.get(i).get(walks.get(i).size() - 1).getKmerAsString());
                        if (gf != null) { Graphs.addGraph(g, gf); }
                    }

                    ConnectivityInspector<CortexVertex, CortexEdge> ci = new ConnectivityInspector<>(g);
                    int index = 0;

                    for (Set<CortexVertex> s : ci.connectedSets()) {
                        String name = String.format("%s:contig%d", parentName.iterator().next(), index);
                        List<CortexVertex> l = TraversalUtils.toWalk(g, s.iterator().next().getKmerAsString(), g.edgeSet().iterator().next().getColor());

                        if (l.size() == 0 && s.size() > 0) {
                            TraversalEngine e = new TraversalEngineFactory()
                                    .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                                    .combinationOperator(OR)
                                    .traversalDirection(BOTH)
                                    .stoppingRule(ContigStopper.class)
                                    .maxBranchLength(2*s.size())
                                    .graph(GRAPH)
                                    .links(LINKS)
                                    .make();

                            for (CortexVertex v : s) {
                                List<CortexVertex> ll = e.walk(v.getKmerAsString());

                                if (ll.size() > l.size()) {
                                    l = ll;
                                }
                            }
                        }

                        if (l.size() > 0) {
                            String cl = TraversalUtils.toContig(l);

                            if (BACKGROUNDS != null) {
                                for (String pn : parentName) {
                                    if (BACKGROUNDS.containsKey(pn)) {
                                        List<SAMRecord> srs = BACKGROUNDS.get(pn).align(cl);

                                        if (srs.size() > 0) {
                                            srs.sort((s1, s2) -> {
                                                int s1length = s1.getAlignmentEnd() - s1.getAlignmentStart();
                                                int nm1 = s1.getIntegerAttribute("NM");

                                                int s2length = s2.getAlignmentEnd() - s2.getAlignmentStart();
                                                int nm2 = s2.getIntegerAttribute("NM");

                                                if (s1length != s2length) {
                                                    return s1length > s2length ? -1 : 1;
                                                }

                                                if (nm1 != nm2) {
                                                    return nm1 < nm2 ? -1 : 1;
                                                }

                                                return 0;
                                            });

                                            SAMRecord sr = srs.get(0);
                                            Interval it = new Interval(sr.getReferenceName(), sr.getAlignmentStart(), sr.getAlignmentEnd(), sr.getReadNegativeStrandFlag(), sr.getCigarString());

                                            String seq = BACKGROUNDS.get(pn).find(it);
                                            targets.put(String.format("%s:%s:%d-%d:%s:%s", pn, it.getContig(), it.getStart(), it.getEnd(), it.isPositiveStrand() ? "+" : "-", it.getName()), seq);
                                        } else {
                                            targets.put(name, cl);
                                        }
                                    }
                                }
                            } else {
                                targets.put(name, cl);
                            }

                            index++;
                        }
                    }

                    if (index == 0) {

                    }
                }

                if (targets.size() > 0) {
                    MosaicAligner ma = new MosaicAligner();
                    List<Triple<String, Pair<Integer, Integer>, String>> lps = ma.align(query, targets);

                    int maxLength = 0;
                    for (Triple<String, Pair<Integer, Integer>, String> lp : lps) {
                        String name = String.format("%s (%d-%d)", lp.getLeft(), lp.getMiddle().getFirst(), lp.getMiddle().getSecond());
                        maxLength = Math.max(maxLength, name.length());
                    }

                    StringBuilder sb = new StringBuilder(StringUtil.repeatCharNTimes(' ', query.length()));
                    for (int i = 0; i <= query.length() - GRAPH.getKmerSize(); i++) {
                        CanonicalKmer ck = new CanonicalKmer(query.substring(i, i + GRAPH.getKmerSize()));

                        if (rois.contains(ck)) {
                            for (int j = i; j <= i + GRAPH.getKmerSize(); j++) {
                                sb.setCharAt(j, '*');
                            }
                        }
                    }

                    for (int i = 0; i < lps.get(0).getRight().length(); i++) {
                        if (lps.get(0).getRight().charAt(i) == '-') {
                            sb.insert(i, sb.charAt(i) == '*' ? '*' : ' ');
                        }
                    }

                    log.info("\n{} {}\n{}", String.format("%" + maxLength + "s", "novel"), sb.toString(), ma);
                }
            }
        }
    }

    private DirectedWeightedPseudograph<CortexVertex, CortexEdge> assembleMiddle(Set<String> parentName, List<CortexVertex> wLeft, List<CortexVertex> wRight) {
        int childColor = GRAPH.getColorForSampleName(ROIS.getSampleName(0));

        for (int i = wLeft.size() - 1; i > 0; i--) {
            Map<Integer, Set<CortexByteKmer>> cbkOut = TraversalUtils.getAllNextKmers(wLeft.get(i).getCortexRecord(), !wLeft.get(i).getKmerAsString().equals(wLeft.get(i).getCanonicalKmer().getKmerAsString()));
            Set<CortexByteKmer> outgoingEdges = new HashSet<>();
            for (int c : GRAPH.getColorsForSampleNames(parentName)) {
                outgoingEdges.addAll(cbkOut.get(c));
            }
            outgoingEdges.removeAll(cbkOut.get(childColor));

            if (outgoingEdges.size() > 0) {
                for (int j = 0; j < wRight.size() - 1; j++) {
                    Map<Integer, Set<CortexByteKmer>> cbkIn = TraversalUtils.getAllPrevKmers(wRight.get(j).getCortexRecord(), !wRight.get(j).getKmerAsString().equals(wRight.get(j).getCanonicalKmer().getKmerAsString()));
                    Set<CortexByteKmer> incomingEdges = new HashSet<>();
                    for (int c : GRAPH.getColorsForSampleNames(parentName)) {
                        incomingEdges.addAll(cbkIn.get(c));
                    }
                    incomingEdges.removeAll(cbkIn.get(childColor));

                    if (incomingEdges.size() > 0) {
                        TraversalEngine ef = new TraversalEngineFactory()
                                .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                                .traversalDirection(FORWARD)
                                .combinationOperator(OR)
                                .stoppingRule(DestinationStopper.class)
                                .maxBranchLength(1000)
                                .graph(GRAPH)
                                .links(LINKS)
                                .make();

                        TraversalEngine er = new TraversalEngineFactory()
                                .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                                .traversalDirection(REVERSE)
                                .combinationOperator(OR)
                                .stoppingRule(DestinationStopper.class)
                                .maxBranchLength(1000)
                                .graph(GRAPH)
                                .links(LINKS)
                                .make();

                        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = ef.dfs(wLeft.get(i-1).getKmerAsString(), wRight.get(j+1).getKmerAsString());
                        if (g == null || g.vertexSet().size() == 0) {
                            g = er.dfs(wRight.get(j+1).getKmerAsString(), wLeft.get(i-1).getKmerAsString());
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

    private Set<CanonicalKmer> loadRois(CortexGraph rg) {
        Set<CanonicalKmer> rois = new HashSet<>();

        for (CortexRecord rr : rg) {
            rois.add(rr.getCanonicalKmer());
        }

        return rois;
    }

    private List<CortexVertex> loadChildWalk(ReferenceSequence seq, CortexGraph graph, List<CortexLinks> links) {
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

        int subcontigStart = regions.get(0).getFirst() - window;
        if (subcontigStart < 0) { subcontigStart = 0; }

        int subcontigStop  = regions.get(regions.size() - 1).getSecond() + window;
        if (subcontigStop >= w.size()) { subcontigStop = w.size() - 1; }

        List<Pair<Integer, Integer>> sections = new ArrayList<>();
        for (int i = 0; i < regions.size() - 1; i++) {
            if (regions.get(i+1).getFirst() - regions.get(i).getSecond() > novelDistanceSplit) {
                sections.add(Pair.create(subcontigStart, regions.get(i).getSecond() + window));

                subcontigStart = regions.get(i+1).getFirst() - window;
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
