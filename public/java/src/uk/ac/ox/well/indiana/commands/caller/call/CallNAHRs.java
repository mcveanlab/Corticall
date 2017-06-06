package uk.ac.ox.well.indiana.commands.caller.call;

import com.google.api.client.util.Joiner;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.jgrapht.traverse.TopologicalOrderIterator;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.stoppingconditions.NahrStopper;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngineFactory;

import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.AND;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

/**
 * Created by kiran on 05/06/2017.
 */
public class CallNAHRs extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="lookups", shortName="l", doc="Lookups")
    public HashMap<String, KmerLookup> LOOKUPS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(PARENTS);

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(GRAPH.getColorForSampleName(CHILD))
                .joiningColors(GRAPH.getColorsForSampleNames(PARENTS))
                .traversalDirection(BOTH)
                .combinationOperator(AND)
                .stopper(NahrStopper.class)
                .graph(GRAPH)
                .rois(ROI)
                .make();

        Map<CortexKmer, Boolean> used = new HashMap<>();
        for (CortexRecord rr : ROI) {
            used.put(rr.getCortexKmer(), false);
        }

        Map<String, IntervalTreeMap<Interval>> candidateLoci = new TreeMap<>();

        for (CortexKmer rk : used.keySet()) {
            if (!used.get(rk)) {
                log.info("{}", rk);

                DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = e.dfs(rk.getKmerAsString());

                TopologicalOrderIterator<CortexVertex, CortexEdge> toi = new TopologicalOrderIterator<>(g);

                Map<String, Integer> refCount = new HashMap<>();
                Map<String, Map<String, Integer>> chrCount = new HashMap<>();
                for (String parent : LOOKUPS.keySet()) {
                    chrCount.put(parent, new HashMap<>());
                }
                int numNovelStretches = 0;
                boolean inNovelStretch = false;

                while (toi.hasNext()) {
                    CortexVertex cv = toi.next();

                    for (String parent : LOOKUPS.keySet()) {
                        Set<Interval> intervals = LOOKUPS.get(parent).findKmer(cv.getSk());

                        if (intervals.size() == 1) {
                            ContainerUtils.increment(refCount, parent);
                            ContainerUtils.increment(chrCount.get(parent), intervals.iterator().next().getContig());

                            Interval newInterval = intervals.iterator().next();

                            if (!candidateLoci.containsKey(newInterval.getContig())) {
                                candidateLoci.put(newInterval.getContig(), new IntervalTreeMap<>());
                            }

                            candidateLoci.get(newInterval.getContig()).put(newInterval, newInterval);
                        }
                    }

                    if (used.containsKey(cv.getCk())) {
                        if (!inNovelStretch) {
                            numNovelStretches++;
                            inNovelStretch = true;
                        }

                        used.put(cv.getCk(), true);
                    } else {
                        inNovelStretch = false;
                    }
                }

                String mostFrequentBackground = mostFrequentBackground(refCount, 3);

                log.info(" -- stats {} {} {} {}", numNovelStretches, refCount, mostFrequentBackground, chrCount);

                if (mostFrequentBackground != null && numNovelStretches > 1 && chrCount.get(mostFrequentBackground).keySet().size() > 1) {
                    toi = new TopologicalOrderIterator<>(g);
                    while (toi.hasNext()) {
                        CortexVertex cv = toi.next();

                        Set<Interval> intervals = LOOKUPS.get(mostFrequentBackground).findKmer(cv.getSk());

                        if (intervals.size() == 1) {
                            log.info(" -- {} {} {} {}", mostFrequentBackground, used.containsKey(cv.getCk()), intervals, recordToString(cv.getCr(), childColor, parentColors));
                        } else {
                            log.info(" -- {} [{}] {} {}", mostFrequentBackground, used.containsKey(cv.getCk()), intervals.size(), recordToString(cv.getCr(), childColor, parentColors));
                        }
                    }
                }

                log.info("");
            }
        }

        for (String contig : candidateLoci.keySet()) {
            for (Interval interval : candidateLoci.get(contig).keySet()) {
                log.info("{}", interval);
            }
        }
    }

    private String mostFrequentBackground(Map<String, Integer> refCount, int thresholdMultiplier) {
        Map.Entry<String, Integer> maxEntry = null;

        for (Map.Entry<String, Integer> e : refCount.entrySet()) {
            if (maxEntry == null || e.getValue() > maxEntry.getValue()) {
                maxEntry = e;
            }
        }

        boolean meetsThreshold = false;

        if (maxEntry != null) {
            for (Map.Entry<String, Integer> e : refCount.entrySet()) {
                if (!maxEntry.getKey().equals(e.getKey()) && maxEntry.getValue() >= thresholdMultiplier*e.getValue()) {
                    meetsThreshold = true;
                }
            }
        }

        return meetsThreshold ? maxEntry.getKey() : null;
    }

    private String recordToString(CortexRecord cr, int childColor, List<Integer> parentColors) {
        List<Object> pieces = new ArrayList<>();
        pieces.add(cr.getKmerAsString());
        pieces.add(cr.getCoverage(childColor));
        parentColors.forEach(c -> pieces.add(cr.getCoverage(c)));
        pieces.add(cr.getEdgesAsString(childColor));
        parentColors.forEach(c -> pieces.add(cr.getEdgesAsString(c)));

        return Joiner.on(' ').join(pieces);
    }
}
