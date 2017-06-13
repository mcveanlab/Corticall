package uk.ac.ox.well.indiana.commands.caller.call;

import com.google.api.client.util.Joiner;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.StringUtil;
import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;
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
import uk.ac.ox.well.indiana.utils.visualizer.GraphVisualizer;

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

        String sk = "ACATGTGGTTCAGGAGAATGGGCTAAAGACAAATGCCGCTGTAAGGA";

        Pair<List<String>, List<Interval>> recon = reconstruct("ref", sk);

        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = buildGraph(recon);

        Map<String, Interval> aggregatedIntervals = aggregateIntervals(mergeIntervals(recon));
        Map<String, Integer> contigIndices = new HashMap<>();

        int index = aggregatedIntervals.size();
        for (String contig : aggregatedIntervals.keySet()) {
            contigIndices.put(contig, index);
            index--;
        }

        List<ReferenceSequence> rseqs = new ArrayList<>();

        for (String contig : aggregatedIntervals.keySet()) {
            Interval it = aggregatedIntervals.get(contig);
            ReferenceSequence rseq = LOOKUPS.get("ref").getReferenceSequence().getSubsequenceAt(it.getContig(), it.getStart(), it.getEnd());
            ReferenceSequence nseq = new ReferenceSequence(rseq.getName(), rseq.getContigIndex(), it.isPositiveStrand() ? rseq.getBaseString().getBytes() : SequenceUtils.reverseComplement(rseq.getBaseString()).getBytes());

            rseqs.add(nseq);
        }

        out.println(Joiner.on('\t').join(Arrays.asList("kmer", "interval", "pos", "contigIndex")));

        for (int i = 0; i < recon.getFirst().size(); i++) {
            String kmer = recon.getFirst().get(i);
            Interval interval = recon.getSecond().get(i);
            int contigIndex = interval == null ? -1 : contigIndices.get(interval.getContig());

            log.info("{} {} {}", kmer, interval, LOOKUPS.get("ref").findKmer(kmer));

            if (contigIndex >= 0) {
                //log.info("{} {}:{}-{},{} {} {}", kmer, interval.getContig(), interval.getStart(), interval.getEnd(), interval.isPositiveStrand() ? "+" : "-", i, contigIndex);
                String intervalString = interval.getContig() + ":" + interval.getStart() + "-" + interval.getEnd() + ":" + (interval.isPositiveStrand() ? "+" : "-");
                out.println(Joiner.on('\t').join(Arrays.asList(kmer, intervalString, i, contigIndex)));
            } else {
                //log.info("{} NA {} {}", kmer, i, contigIndex);
                out.println(Joiner.on('\t').join(Arrays.asList(kmer, "NA", i, contigIndex)));
            }
        }

        GraphVisualizer gv = new GraphVisualizer(9000);
        gv.display(g, rseqs, "nahr1");
    }

    private DirectedWeightedPseudograph<CortexVertex, CortexEdge> buildGraph(Pair<List<String>, List<Interval>> recon) {
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = new DirectedWeightedPseudograph<>(CortexEdge.class);

        for (int i = 0; i < recon.getFirst().size() - 1; i++) {
            String s0 = recon.getFirst().get(i);
            CortexRecord c0 = GRAPH.findRecord(new CortexKmer(s0));
            Interval i0 = recon.getSecond().get(i);
            CortexVertex v0 = new CortexVertex(s0, c0, i0);

            String s1 = recon.getFirst().get(i + 1);
            CortexRecord c1 = GRAPH.findRecord(new CortexKmer(s1));
            Interval i1 = recon.getSecond().get(i + 1);
            CortexVertex v1 = new CortexVertex(s1, c1, i1);

            g.addVertex(v0);
            g.addVertex(v1);

            Map<Integer, Set<String>> aks = TraversalEngine.getAllNextKmers(c0, !s0.equals(c0.getKmerAsString()));
            Set<String> pnk = aks.get(GRAPH.getColorForSampleName("ref"));
            Set<String> cnk = aks.get(GRAPH.getColorForSampleName(CHILD));

            if (pnk.contains(s1)) {
                g.addEdge(v0, v1, new CortexEdge(GRAPH.getColorForSampleName("ref"), 1.0));
            }

            if (cnk.contains(s1)) {
                g.addEdge(v0, v1, new CortexEdge(GRAPH.getColorForSampleName(CHILD), 1.0));
            }
        }

        return g;
    }

    @NotNull
    private Map<String, Interval> aggregateIntervals(Set<Interval> mergedIntervals) {
        Map<String, Interval> aggregatedIntervals = new TreeMap<>();
        for (Interval it : mergedIntervals) {
            if (!aggregatedIntervals.containsKey(it.getContig())) {
                aggregatedIntervals.put(it.getContig(), it);
            } else {
                Interval locus = aggregatedIntervals.get(it.getContig());

                int start = locus.getStart() < it.getStart() ? locus.getStart() : it.getStart();
                int end   = locus.getEnd()   > it.getEnd()   ? locus.getEnd()   : it.getEnd();

                Interval newlocus = new Interval(locus.getContig(), start, end, locus.isNegativeStrand(), null);

                aggregatedIntervals.put(locus.getContig(), newlocus);
            }
        }
        return aggregatedIntervals;
    }

    @NotNull
    private Set<Interval> mergeIntervals(Pair<List<String>, List<Interval>> recon) {
        Set<Interval> mergedIntervals = new TreeSet<>();

        Interval locus = null;
        for (Interval it : recon.getSecond()) {
            if (locus == null) {
                locus = it;
            } else {
                if (it != null && locus.getContig().equals(it.getContig()) && locus.isPositiveStrand() == it.isPositiveStrand()) {
                    int start = locus.getStart() < it.getStart() ? locus.getStart() : it.getStart();
                    int end   = locus.getEnd()   > it.getEnd()   ? locus.getEnd()   : it.getEnd();

                    locus = new Interval(locus.getContig(), start, end, locus.isNegativeStrand(), null);
                } else {
                    mergedIntervals.add(locus);

                    locus = it;
                }
            }
        }

        if (locus != null) {
            mergedIntervals.add(locus);
        }

        return mergedIntervals;
    }

    private Pair<List<String>, List<Interval>> reconstruct(String background, String sk) {
        Pair<List<String>, List<Interval>> rev = reconstruct(background, sk, false, 5000);
        Pair<List<String>, List<Interval>> fwd = reconstruct(background, sk, true, 5000);

        List<String> allKmers = new ArrayList<>();
        List<Interval> allLoci = new ArrayList<>();

        allKmers.addAll(rev.getFirst());
        allKmers.add(sk);
        allKmers.addAll(fwd.getFirst());

        allLoci.addAll(rev.getSecond());
        allLoci.add(null);
        allLoci.addAll(fwd.getSecond());

        return new Pair<>(allKmers, allLoci);
    }

    private Pair<List<String>, List<Interval>> reconstruct(String background, String sk, boolean goForward, int limit) {
        List<String> vertices = new ArrayList<>();
        List<Interval> loci = new ArrayList<>();

        Interval ci = null;
        int distanceFromNovel = 0;
        boolean onRef = false;
        boolean positiveStrand = false;
        boolean keepGoing = true;

        while (distanceFromNovel < limit && keepGoing) {
            if (ci == null) {
                CortexRecord cr = GRAPH.findRecord(new CortexKmer(sk));
                Map<Integer, Set<String>> aks = goForward ? TraversalEngine.getAllNextKmers(cr, !sk.equals(cr.getKmerAsString())) : TraversalEngine.getAllPrevKmers(cr, !sk.equals(cr.getKmerAsString()));
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
                    //log.info("{} {} {} {} {}", vertices.get(vertices.size() - 1), loci.get(loci.size() - 1), vertices.get(0), loci.get(0), recordToString(cr, GRAPH.getColorForSampleName(CHILD), GRAPH.getColorsForSampleNames(PARENTS)));

                    if (ci != null) {
                        onRef = true;
                        positiveStrand = ci.isPositiveStrand();
                    }
                } else {
                    //log.info("lastNovel onRef={} achi.size()={}", onRef, achi.size());
                    keepGoing = false;
                }
            } else {
                do {
                    Interval aci;

                    if (positiveStrand) {
                        if (goForward) {
                            aci = new Interval(ci.getContig(), ci.getStart() + 1, ci.getEnd() + 1, ci.isNegativeStrand(), null);
                        } else {
                            aci = new Interval(ci.getContig(), ci.getStart() - 1, ci.getEnd() - 1, ci.isNegativeStrand(), null);
                        }
                    } else {
                        if (goForward) {
                            aci = new Interval(ci.getContig(), ci.getStart() - 1, ci.getEnd() - 1, ci.isNegativeStrand(), null);
                        } else {
                            aci = new Interval(ci.getContig(), ci.getStart() + 1, ci.getEnd() + 1, ci.isNegativeStrand(), null);
                        }
                    }

                    String aref = LOOKUPS.get(background).findKmer(aci);

                    CortexRecord cr = GRAPH.findRecord(new CortexKmer(sk));
                    Map<Integer, Set<String>> aks = goForward ? TraversalEngine.getAllNextKmers(cr, !sk.equals(cr.getKmerAsString())) : TraversalEngine.getAllPrevKmers(cr, !sk.equals(cr.getKmerAsString()));
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
                        //log.info("{} {} {} {}", vertices.get(vertices.size() - 1), loci.get(loci.size() - 1), vertices.get(0), loci.get(0));
                        //log.info("{} {} {} {} {}", vertices.get(vertices.size() - 1), loci.get(loci.size() - 1), vertices.get(0), loci.get(0), recordToString(cr, GRAPH.getColorForSampleName(CHILD), GRAPH.getColorsForSampleNames(PARENTS)));
                    } else {
                        if (achi.size() == 1) {
                            distanceFromNovel = 0;
                            sk = achi.iterator().next();
                            Set<Interval> acis = LOOKUPS.get(background).findKmer(sk);
                            ci = acis.size() == 1 ? acis.iterator().next() : null;

                            if ((ci != null && !aci.getContig().equals(ci.getContig())) || acis.size() > 1) {
                                keepGoing = false;
                            } else if (goForward) {
                                vertices.add(sk);
                                loci.add(ci);
                            } else {
                                vertices.add(0, sk);
                                loci.add(0, ci);
                            }

                            //log.info("{} {} {} {}", vertices.get(vertices.size() - 1), loci.get(loci.size() - 1), vertices.get(0), loci.get(0));
                            //log.info("{} {} {} {} {}", vertices.get(vertices.size() - 1), loci.get(loci.size() - 1), vertices.get(0), loci.get(0), recordToString(cr, GRAPH.getColorForSampleName(CHILD), GRAPH.getColorsForSampleNames(PARENTS)));
                        } else {
                            //log.info("firstNovel onRef={} achi.size()={}", onRef, achi.size());
                            keepGoing = false;
                        }

                        onRef = false;
                    }
                } while (onRef);
            }
        }

        //log.info("  {}", vertices.size());
        //log.info("  {}", loci.size());
        //log.info("");

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
