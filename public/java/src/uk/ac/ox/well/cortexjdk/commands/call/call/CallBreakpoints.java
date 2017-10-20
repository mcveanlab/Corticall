package uk.ac.ox.well.cortexjdk.commands.call.call;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Interval;
import org.jgrapht.DirectedGraph;
import org.jgrapht.GraphPath;
import org.jgrapht.Graphs;
import org.jgrapht.alg.shortestpath.DijkstraShortestPath;
import org.jgrapht.alg.util.Pair;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.jgrapht.traverse.DepthFirstIterator;
import org.jgrapht.traverse.TopologicalOrderIterator;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.alignment.pairwise.BwaAligner;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.caller.Bubble;
import uk.ac.ox.well.cortexjdk.utils.caller.BubbleCaller;
import uk.ac.ox.well.cortexjdk.utils.caller.BubbleCallerFactory;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.ConnectivityAnnotations;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.BubbleClosingStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.BubbleOpeningStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.NovelContinuationStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

import static org.jgrapht.Graphs.successorListOf;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.AND;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.REVERSE;

/**
 * Created by kiran on 30/08/2017.
 */
public class CallBreakpoints extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="links", shortName="l", doc="Links")
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName="references", shortName="R", doc="References")
    public HashMap<String, KmerLookup> REFERENCES;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="mq", shortName="mq", doc="Mapping quality minimum threshold")
    public Integer MQ_THRESHOLD = 10;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CortexKmer, Boolean> seen = new HashMap<>();
        for (CortexRecord cr : ROI) {
            seen.put(cr.getCortexKmer(), false);
        }

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(GRAPH.getColorForSampleName(ROI.getSampleName(0)))
                .joiningColors(GRAPH.getColorsForSampleNames(REFERENCES.keySet()))
                .combinationOperator(OR)
                .stoppingRule(NovelContinuationStopper.class)
                .graph(GRAPH)
                .links(LINKS)
                .references(REFERENCES.values())
                .rois(ROI)
                .make();

        for (CortexKmer ck : seen.keySet()) {
            if (!seen.get(ck)) {
                List<CortexVertex> w = e.walk(ck.getKmerAsString());

                int wOldSize = w.size();

                boolean extended = false;
                do {
                    extended = false;
                    List<List<CortexVertex>> extFwd = new ArrayList<>();

                    Set<CortexVertex> nvs = e.getNextVertices(w.get(w.size() - 1).getBk());
                    for (CortexVertex cv : nvs) {
                        List<CortexVertex> wn = e.walk(cv.getSk(), true);
                        wn.add(0, cv);

                        boolean hasNovels = false;

                        for (CortexVertex v : wn) {
                            if (seen.containsKey(v.getCk()) && !seen.get(v.getCk())) {
                                hasNovels = true;
                                break;
                            }
                        }

                        if (hasNovels) {
                            extFwd.add(wn);
                        }
                    }

                    if (extFwd.size() == 1) {
                        w.addAll(extFwd.get(0));
                        extended = true;

                        for (CortexVertex v : extFwd.get(0)) {
                            if (seen.containsKey(v.getCk())) {
                                seen.put(v.getCk(), true);
                            }
                        }
                    }
                } while (extended);

                do {
                    extended = false;
                    List<List<CortexVertex>> extRev = new ArrayList<>();

                    Set<CortexVertex> pvs = e.getPrevVertices(w.get(0).getBk());
                    for (CortexVertex cv : pvs) {
                        List<CortexVertex> wp = e.walk(cv.getSk(), false);
                        wp.add(cv);

                        boolean hasNovels = false;

                        for (CortexVertex v : wp) {
                            if (seen.containsKey(v.getCk()) && !seen.get(v.getCk())) {
                                hasNovels = true;
                                break;
                            }
                        }

                        if (hasNovels) {
                            extRev.add(wp);
                        }
                    }

                    if (extRev.size() == 1) {
                        w.addAll(0, extRev.get(0));
                        extended = true;

                        for (CortexVertex v : extRev.get(0)) {
                            if (seen.containsKey(v.getCk())) {
                                seen.put(v.getCk(), true);
                            }
                        }
                    }
                } while (extended);

                int wNewSize = w.size();

                log.info("{} {} {}", ck, wOldSize, wNewSize);

                /*
                log.info("{} {} {} {}", ck, w.size(), bor == null ? null : bor.vertexSet().size(), bof == null ? null : bof.vertexSet().size());

                e.getConfiguration().setTraversalDirection(REVERSE);
                DirectedWeightedPseudograph<CortexVertex, CortexEdge> bor = e.dfs(w.get(0).getSk());

                Set<CortexVertex> nvs = e.getNextVertices(w.get(w.size() - 1).getBk());

                e.getConfiguration().setTraversalDirection(FORWARD);
                for (CortexVertex cv : nvs) {
                    DirectedWeightedPseudograph<CortexVertex, CortexEdge> bof = e.dfs(cv.getSk());
                }

                log.info("{} {} {} {}", ck, w.size(), bor == null ? null : bor.vertexSet().size(), bof == null ? null : bof.vertexSet().size());

                if (ck.equals(new CortexKmer("TAGAACCACGTTTTTCGGATACACTTGGTGGCCAGTGTACTAACAAA"))) {
                    for (int i = w.size() - 10; i < w.size(); i++) {
                        log.info("  {} {}", i, w.get(i));
                    }

                    TopologicalOrderIterator<CortexVertex, CortexEdge> toi = new TopologicalOrderIterator<>(bof);
                    while (toi.hasNext()) {
                        CortexVertex cn = toi.next();
                        log.info("  {} {} {}", bof.inDegreeOf(cn), bof.outDegreeOf(cn), cn);
                    }

                    log.info("");
                }
                */

                /*
                for (String parent : REFERENCES.keySet()) {
                    log.info("{} {} {}", ck, parent, w.size());

                    List<CortexVertex> p = closeBubbles(w, parent, seen);

//                    log.info("  w={}", TraversalEngine.toContig(w));
//                    log.info("  p={}", TraversalEngine.toContig(p));
//
//                    List<SAMRecord> srs = REFERENCES.get(parent).getAligner().align(TraversalEngine.toContig(p));
//
//                    for (SAMRecord sr : srs) {
//                        if (sr.getMappingQuality() > 0) {
//                            log.info("  {}", sr.getSAMString().trim());
//                        }
//                    }
//
//                    log.info("");
                }
                */

                //for (CortexVertex v : bo.vertexSet()) {
                for (CortexVertex v : w) {
                    if (seen.containsKey(v.getCk())) {
                        seen.put(v.getCk(), true);
                    }
                }
            }
        }
    }

    private List<CortexVertex> closeBubbles(List<CortexVertex> w, String parent, Map<CortexKmer, Boolean> seen) {
        //Map<Integer, DirectedWeightedPseudograph<CortexVertex, CortexEdge>> bubbles = new TreeMap<>();

        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = new DirectedWeightedPseudograph<>(CortexEdge.class);
        Map<CortexVertex, Integer> indices = new HashMap<>();
        indices.put(w.get(0), 0);

        for (int i = 1; i < w.size(); i++) {
            CortexVertex v0 = w.get(i - 1);
            CortexVertex v1 = w.get(i);

            g.addVertex(v0);
            g.addVertex(v1);
            g.addEdge(v0, v1, new CortexEdge());

            indices.put(w.get(i), i);
        }

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(GRAPH.getColorForSampleName(parent))
                .joiningColors(GRAPH.getColorForSampleName(ROI.getSampleName(0)))
                .recruitmentColors(GRAPH.getColorForSampleName(REFERENCES.get(parent).getSources().iterator().next()))
                .traversalDirection(FORWARD)
                .combinationOperator(OR)
                //.connectAllNeighbors(true)
                .stoppingRule(BubbleClosingStopper.class)
                .previousTraversal(g)
                .graph(GRAPH)
                .links(LINKS)
                .references(REFERENCES.values())
                .make();

        for (int i = 0; i < w.size() - 1; i++) {
            CortexVertex vi = w.get(i);

            if (seen.containsKey(vi.getCk())) {
                log.info("* {} {} novel", i, vi);

                Set<CortexVertex> sources = new LinkedHashSet<>();

                int lowerLimit = i - vi.getCr().getKmerSize() >= 0 ? i - vi.getCr().getKmerSize() : 0;
                for (int j = i - 1; j >= lowerLimit; j--) {
                    CortexVertex vj = w.get(j);
                    CortexVertex vk = w.get(j+1);

                    if (!seen.containsKey(vj.getCk())) {
                        Set<CortexVertex> nvs = e.getNextVertices(new CortexByteKmer(vj.getSk()));

                        for (CortexVertex cv : nvs) {
                            if (!cv.equals(vk)) {
                                sources.add(cv);
                            }
                        }
                    }
                }

                for (CortexVertex source : sources) {
                    DirectedWeightedPseudograph<CortexVertex, CortexEdge> b = e.dfs(source.getSk());

                    if (b != null) {
                        CortexVertex sink = b.vertexSet().iterator().next();
                        for (CortexVertex v : b.vertexSet()) {
                            log.info("  {} {}", v.getCk(), indices.get(v));
                        }

                        log.info("Hello!");
                    }
                }
            }
        }

        /*
        Map<CortexVertex, Integer> sources = new HashMap<>();
        Map<CortexVertex, Integer> sinks = new HashMap<>();

        for (int i = 0; i < w.size() - 1; i++) {
            CortexVertex v0 = w.get(i);
            CortexVertex v1 = w.get(i + 1);

            Set<CortexVertex> nvs = e.getNextVertices(new CortexByteKmer(v0.getSk()));

            for (CortexVertex cv : nvs) {
                if (!cv.equals(v1)) {
                    sources.put(cv, i);
                }
            }

            Set<CortexVertex> pvs = e.getPrevVertices(new CortexByteKmer(v1.getSk()));

            for (CortexVertex cv : pvs) {
                if (!cv.equals(v0)) {
                    sinks.put(cv, i + 1);
                }
            }
        }

        List<CortexVertex> wp = new ArrayList<>();
        for (int i = 0; i < w.size(); i++) {
            if (bubbles.containsKey(i)) {
                DirectedWeightedPseudograph<CortexVertex, CortexEdge> b = bubbles.get(i);

                CortexVertex source = w.get(i);
                CortexVertex sink = null;
                int sinkIndex = i;

                for (CortexVertex v : b.vertexSet()) {
                    if (indices.containsKey(v) && indices.get(v) > sinkIndex) {
                        sink = v;
                        sinkIndex = indices.get(v);
                    }
                }

                if (sink != null) {
                    GraphPath<CortexVertex, CortexEdge> gp = DijkstraShortestPath.findPathBetween(b, source, sink);
                    List<CortexVertex> vp = gp.getVertexList();

                    for (CortexVertex v : vp) {
                        wp.add(v);
                    }

                    int j = w.indexOf(sink);
                    if (j == -1) {
                        throw new CortexJDKException("Sink vertex for bubble is missing from contig graph.");
                    }

                    i = j;
                } else {
                    wp.add(w.get(i));
                }
            } else {
                wp.add(w.get(i));
            }
        }

        return wp;
        */

        return null;
    }

    private SAMRecord chooseBestAlignment(String contig, int mqThreshold) {
        SAMRecord bestRecord = null;
        int bestRecordScore = Integer.MAX_VALUE;

        for (KmerLookup kl : REFERENCES.values()) {
            List<SAMRecord> recs = kl.getAligner().align(contig);
            List<SAMRecord> filteredRecs = new ArrayList<>();

            for (SAMRecord rec : recs) {
                if (rec.getMappingQuality() >= mqThreshold) {
                    filteredRecs.add(rec);
                }
            }

            if (filteredRecs.size() == 1) {
                SAMRecord filteredRec = filteredRecs.get(0);
                int filteredRecScore = scoreAlignment(filteredRec);

                if (filteredRecScore < bestRecordScore) {
                    bestRecord = filteredRec;
                    bestRecordScore = filteredRecScore;
                }
            }
        }

        return bestRecord;
    }

    private int scoreAlignment(SAMRecord record) {
        int basesChanged = record.getIntegerAttribute("NM");

        for (CigarElement ce : record.getCigar()) {
            if (ce.getOperator().equals(CigarOperator.SOFT_CLIP) ||
                ce.getOperator().equals(CigarOperator.INSERTION) ||
                ce.getOperator().equals(CigarOperator.DELETION)) {
                basesChanged += ce.getLength();
            }
        }

        return basesChanged;
    }
}
