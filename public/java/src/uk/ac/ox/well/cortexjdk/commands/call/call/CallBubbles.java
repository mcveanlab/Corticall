package uk.ac.ox.well.cortexjdk.commands.call.call;

import com.google.common.base.Joiner;
import htsjdk.samtools.util.Interval;
import org.jgrapht.GraphPath;
import org.jgrapht.alg.shortestpath.DijkstraShortestPath;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.jgrapht.graph.EdgeReversedGraph;
import org.jgrapht.traverse.DepthFirstIterator;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.caller.Bubble;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraphWriter;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.BubbleClosingStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.BubbleOpeningStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.AND;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;

/**
 * Created by kiran on 29/08/2017.
 */
public class CallBubbles extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="links", shortName="l", doc="Links")
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="references", shortName="R", doc="References")
    public HashMap<String, KmerLookup> REFERENCES;

    @Output
    public PrintStream out;

    @Output(fullName="rout", shortName="ro", doc="Remaining ROI out")
    public File rout;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(new ArrayList<>(REFERENCES.keySet()));

        TraversalEngine eOpen = new TraversalEngineFactory()
                .traversalColor(childColor)
                .joiningColors(parentColors)
                .traversalDirection(BOTH)
                .combinationOperator(AND)
                .stoppingRule(BubbleOpeningStopper.class)
                .graph(GRAPH)
                .links(LINKS)
                .rois(ROI)
                .make();

        List<TraversalEngine> eCloses = new ArrayList<>();
        for (int pc : parentColors) {
            TraversalEngine eClose = new TraversalEngineFactory()
                    .traversalColor(parentColors.get(pc))
                    .joiningColors(childColor)
                    .traversalDirection(FORWARD)
                    .combinationOperator(OR)
                    .stoppingRule(BubbleClosingStopper.class)
                    .graph(GRAPH)
                    .links(LINKS)
                    .rois(ROI)
                    .make();

            eCloses.add(eClose);
        }

        out.println(Joiner.on("\t").join("contig", "start", "type", "ref", "alt", "flank5p", "flank3p", "nkCount", "nk", "nks"));

        Map<CortexKmer, Boolean> used = loadRois();
        int numBubbles = 0;
        int numNovelKmersInVariants = 0;

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing novel kmers")
                .message("records processed")
                .maxRecord(used.size())
                .make(log);

        for (CortexKmer rk : used.keySet()) {
            if (!used.get(rk)) {
                DirectedWeightedPseudograph<CortexVertex, CortexEdge> gc = eOpen.dfs(rk.getKmerAsString());

                DepthFirstIterator<CortexVertex, CortexEdge> dFwd = null;
                DepthFirstIterator<CortexVertex, CortexEdge> dRev = null;

                for (CortexVertex cv : gc.vertexSet()) {
                    if (cv.getCk().equals(rk)) {
                        dFwd = new DepthFirstIterator<>(gc, cv);
                        dRev = new DepthFirstIterator<>(new EdgeReversedGraph<>(gc), cv);
                    }
                }

                Set<CortexVertex> sources = getCandidates(dRev);
                Set<CortexVertex> sinks = getCandidates(dFwd);

                Set<Bubble> bubbles = new HashSet<>();
                for (CortexVertex so : sources) {
                    for (CortexVertex si : sinks) {
                        DirectedWeightedPseudograph<CortexVertex, CortexEdge> s = new DirectedWeightedPseudograph<>(CortexEdge.class);
                        s.addVertex(si);

                        for (int pc : parentColors) {
                            Set<Interval> soIntervals = REFERENCES.get(GRAPH.getSampleName(pc)).findKmer(so.getSk());
                            Set<Interval> siIntervals = REFERENCES.get(GRAPH.getSampleName(pc)).findKmer(si.getSk());

                            if (so.getCr().getCoverage(pc) > 0 && si.getCr().getCoverage(pc) > 0 && soIntervals.size() == 1 && siIntervals.size() == 1) {
                                Interval soInterval = soIntervals.iterator().next();
                                Interval siInterval = siIntervals.iterator().next();

                                if (soInterval.getContig().equals(siInterval.getContig())) {
                                    eCloses.get(pc).getConfiguration().setPreviousTraversal(s);

                                    DirectedWeightedPseudograph<CortexVertex, CortexEdge> gp = eCloses.get(pc).dfs(so.getSk());

                                    if (gp != null) {
                                        GraphPath<CortexVertex, CortexEdge> dgc = DijkstraShortestPath.findPathBetween(gc, so, si);
                                        GraphPath<CortexVertex, CortexEdge> dgp = DijkstraShortestPath.findPathBetween(gp, so, si);

                                        Set<CortexKmer> novelKmersInBubble = new HashSet<>();
                                        for (CortexVertex cv : dgc.getVertexList()) {
                                            if (used.containsKey(cv.getCk())) {
                                                novelKmersInBubble.add(cv.getCk());
                                            }
                                        }

                                        Bubble b = new Bubble(dgp, dgc, REFERENCES.get(GRAPH.getSampleName(pc)), novelKmersInBubble);
                                        bubbles.add(b);
                                    }
                                }
                            }
                        }
                    }
                }

                for (Bubble b : bubbles) {
                    for (CortexKmer ck : b.getNovelKmers()) {
                        used.put(ck, true);
                        numNovelKmersInVariants++;
                    }

                    out.println(Joiner.on("\t").join(b.getLocus().getContig(), b.getLocus().getStart(), b.getType(), "'" + b.getRefAllele() + "'", "'" + b.getAltAllele() + "'", b.getFlank5p(), b.getFlank3p(), b.getNovelKmers().size(), rk, Joiner.on(",").join(b.getNovelKmers())));

                    numBubbles++;
                }
            }

            pm.update();
        }

        log.info("Found {} bubbles.  Used {}/{} novel kmers", numBubbles, numNovelKmersInVariants, used.size());

        CortexGraphWriter cgw = new CortexGraphWriter(rout);
        cgw.setHeader(ROI.getHeader());

        for (CortexRecord rr : ROI) {
            if (!used.get(rr.getCortexKmer())) {
                cgw.addRecord(rr);
            }
        }

        cgw.close();
    }

    private Set<CortexVertex> getCandidates(DepthFirstIterator<CortexVertex, CortexEdge> d) {
        Set<CortexVertex> candidates = new HashSet<>();
        if (d != null) {
            while (d.hasNext()) {
                CortexVertex cv = d.next();
                Set<Interval> allIntervals = new HashSet<>();
                for (KmerLookup kl : REFERENCES.values()) {
                    Set<Interval> intervals = kl.findKmer(cv.getSk());
                    if (intervals != null) {
                        allIntervals.addAll(intervals);
                    }
                }

                if (allIntervals.size() == 1) {
                    candidates.add(cv);
                }
            }
        }
        return candidates;
    }

    private Map<CortexKmer, Boolean> loadRois() {
        Map<CortexKmer, Boolean> rrs = new HashMap<>();
        for (CortexRecord rr : ROI) {
            rrs.put(rr.getCortexKmer(), false);
        }
        return rrs;
    }
}
