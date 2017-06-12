package uk.ac.ox.well.indiana.commands.caller.call;

import com.google.api.client.util.Joiner;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.StringUtil;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.jgrapht.traverse.TopologicalOrderIterator;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;
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

        /*
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

        candidateLoci = mergeIntervals(candidateLoci, 1);
        candidateLoci = mergeIntervals(candidateLoci, 10);
        candidateLoci = mergeIntervals(candidateLoci, 100);
        candidateLoci = mergeIntervals(candidateLoci, 1000);
        candidateLoci = mergeIntervals(candidateLoci, 10000);
        candidateLoci = mergeIntervals(candidateLoci, 20000);

        for (String contig : candidateLoci.keySet()) {
            for (Interval interval : candidateLoci.get(contig).keySet()) {
                log.info("{}", interval);
            }
        }
        */

        Interval it = new Interval("Pf3D7_01_v3", 29500, 29546);

        String qk = LOOKUPS.get("ref").findKmer(it);
        String rk = SequenceUtils.reverseComplement(qk);
        Set<Interval> iqs = LOOKUPS.get("ref").findKmer(qk);
        Set<Interval> irs = LOOKUPS.get("ref").findKmer(rk);

        reconstruct("ref", qk);

        //reconstruct("ref", new Interval("Pf3D7_01_v3", 29500, 29600), candidateLoci, used);
    }

    private void reconstruct(String background, String sk) {
        Pair<List<String>, List<Interval>> fwd = reconstruct("ref", sk, true, 5000);
        //Pair<List<String>, List<Interval>> rev = reconstruct("ref", sk, false);
    }

    private Pair<List<String>, List<Interval>> reconstruct(String background, String sk, boolean goForward, int limit) {
        List<String> vertices = new ArrayList<>();
        List<Interval> loci = new ArrayList<>();

        Set<Interval> intervals = LOOKUPS.get(background).findKmer(sk);
        Interval ci = null;
        if (intervals.size() == 1) { ci = intervals.iterator().next(); }

        vertices.add(sk);
        loci.add(ci);
        log.info("{} {}", vertices.get(vertices.size() - 1), loci.get(loci.size() - 1));

        boolean onRef = true;
        int distanceFromNovel = 0;
        boolean positiveStrand = true;

        while (distanceFromNovel < limit) {
            if (ci != null) {
                do {
                    Interval aci = positiveStrand ? new Interval(ci.getContig(), ci.getStart() + 1, ci.getEnd() + 1, ci.isNegativeStrand(), null) : new Interval(ci.getContig(), ci.getStart() - 1, ci.getEnd() - 1, ci.isNegativeStrand(), null);
                    String aref = LOOKUPS.get(background).findKmer(aci);

                    CortexRecord cr = GRAPH.findRecord(new CortexKmer(sk));
                    Map<Integer, Set<String>> aks = goForward ? TraversalEngine.getAllNextKmers(cr, !sk.equals(cr.getKmerAsString())) : TraversalEngine.getAllPrevKmers(cr, !sk.equals(cr.getKmerAsString()));
                    //Map<Integer, Set<String>> aks = TraversalEngine.getAllNextKmers(cr, !sk.equals(cr.getKmerAsString()));
                    Set<String> achi = aks.get(GRAPH.getColorForSampleName(CHILD));

                    if (achi.contains(aref)) {
                        distanceFromNovel++;
                        sk = aref;
                        ci = aci;

                        if (goForward) {
                            vertices.add(sk);
                            loci.add(ci);
                        } else {
                            vertices.add(0, sk);
                            loci.add(0, ci);
                        }
                        log.info("{} {}", vertices.get(vertices.size() - 1), loci.get(loci.size() - 1));
                    } else {
                        if (achi.size() == 1) {
                            distanceFromNovel = 0;
                            sk = achi.iterator().next();
                            ci = null;

                            if (goForward) {
                                vertices.add(sk);
                                loci.add(ci);
                            } else {
                                vertices.add(0, sk);
                                loci.add(0, ci);
                            }
                            log.info("{} {}", vertices.get(vertices.size() - 1), loci.get(loci.size() - 1));
                        } else {
                            log.info("firstNovel onRef={} achi.size()={}", onRef, achi.size());
                            break;
                        }

                        onRef = false;
                    }
                } while (onRef);
            } else {
                CortexRecord cr = GRAPH.findRecord(new CortexKmer(sk));
                Map<Integer, Set<String>> aks = goForward ? TraversalEngine.getAllNextKmers(cr, !sk.equals(cr.getKmerAsString())) : TraversalEngine.getAllPrevKmers(cr, !sk.equals(cr.getKmerAsString()));
                //Map<Integer, Set<String>> aks = TraversalEngine.getAllNextKmers(cr, !sk.equals(cr.getKmerAsString()));
                Set<String> apar = aks.get(GRAPH.getColorForSampleName(background));
                Set<String> achi = aks.get(GRAPH.getColorForSampleName(CHILD));

                if (achi.size() == 1) {
                    distanceFromNovel = 0;
                    sk = achi.iterator().next();

                    Set<Interval> acis = LOOKUPS.get(background).findKmer(sk);
                    ci = acis.size() == 1 ? acis.iterator().next() : null;

                    if (goForward) {
                        vertices.add(sk);
                        loci.add(ci);
                    } else {
                        vertices.add(0, sk);
                        loci.add(0, ci);
                    }
                    log.info("{} {}", vertices.get(vertices.size() - 1), loci.get(loci.size() - 1));

                    if (ci != null) {
                        onRef = true;
                        //goForward = ci.isPositiveStrand();
                        positiveStrand = ci.isPositiveStrand();
                    }
                } else {
                    log.info("lastNovel onRef={} achi.size()={}", onRef, achi.size());
                    break;
                }
            }
        }

        log.info("  {}", vertices.size());
        log.info("  {}", loci.size());
        log.info("");

        return new Pair<>(vertices, loci);
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
