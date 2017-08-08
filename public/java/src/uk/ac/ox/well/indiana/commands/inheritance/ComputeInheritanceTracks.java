package uk.ac.ox.well.indiana.commands.inheritance;

import htsjdk.samtools.util.Interval;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexBinaryKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.stoppingconditions.BubbleClosingStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.BubbleOpeningStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.DestinationStopper;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngineFactory;

import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.AND;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;

/**
 * Created by kiran on 08/08/2017.
 */
public class ComputeInheritanceTracks extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="parent", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="drafts", shortName="d", doc="Parental drafts")
    public HashMap<String, KmerLookup> DRAFTS;

    @Argument(fullName="reference", shortName="R", doc="Canonical reference")
    public KmerLookup REFERENCE;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Set<Integer> parentColors = new TreeSet<>(GRAPH.getColorsForSampleNames(PARENTS));
        Set<Integer> draftColors = new TreeSet<>(GRAPH.getColorsForSampleNames(new ArrayList<>(DRAFTS.keySet())));
        int refColor = GRAPH.getColorForSampleName("ref");
        Set<Integer> childColors = getChildColors(parentColors, draftColors, refColor);

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing records")
                .message("records processed")
                .maxRecord(GRAPH.getNumRecords())
                .make(log);

        Set<CortexBinaryKmer> seen = new HashSet<>();

        for (CortexRecord cr : GRAPH) {
            if (!seen.contains(cr.getCortexBinaryKmer()) &&
                isSharedWithOnlyOneParent(cr, parentColors, draftColors) &&
                isSharedWithSomeChildren(cr, parentColors, draftColors, refColor) &&
                isSinglyConnected(cr) &&
                hasUniqueCoordinates(cr, draftColors) &&
                canWalkToCanonicalReference(cr, childColors, refColor)) {

                List<Interval> intervals = getCanonicalReferenceCoordinates(cr, childColors, refColor);
                int draftColor = getDraftColor(cr, draftColors);

                if (intervals != null) {
                    int bubbleColor = getBubbleColor(draftColors, draftColor);

                    Set<String> childAllele = new HashSet<>();
                    Set<String> draftAllele = new HashSet<>();
                    Set<Interval> locus = new HashSet<>();
                    Set<Integer> colors = new HashSet<>();

                    for (int cc : childColors) {
                        if (cr.getCoverage(cc) > 0) {
                            Pair<String, String> bubble = callSimpleBubble(cr, cc, bubbleColor);

                            if (bubble != null) {
                                Interval coords = getBubbleCanonicalCoordinates(bubble.getFirst(), cc, refColor);
                                Pair<String, String> alleles = trimToAlleles(bubble);

                                //log.info("  - {} {} {} {}", cr.getKmerAsString(), alleles.getFirst(), alleles.getSecond(), coords);

                                childAllele.add(alleles.getFirst());
                                draftAllele.add(alleles.getSecond());
                                locus.add(coords);
                                colors.add(cc);

                                for (int i = 0; i <= bubble.getFirst().length() - GRAPH.getKmerSize(); i++) {
                                    String sk = bubble.getFirst().substring(i, i + GRAPH.getKmerSize());
                                    seen.add(new CortexBinaryKmer(CortexRecord.encodeBinaryKmer(sk.getBytes())));
                                }
                            }
                        }
                    }

                    if (childAllele.size() == 1) {
                        log.info("  {} {} {} {}", childAllele, draftAllele, locus, colors);
                    }
                }
            }

            pm.update("records processed, " + seen.size() + " kmers from variants");
        }
    }

    private int getBubbleColor(Set<Integer> draftColors, int draftColor) {
        int bubbleColor = draftColor;
        for (int dc : draftColors) {
            if (dc != draftColor) {
                bubbleColor = dc;
            }
        }
        return bubbleColor;
    }

    private Set<Integer> getChildColors(Set<Integer> parentColors, Set<Integer> draftColors, int refColor) {
        Set<Integer> childColors = new TreeSet<>();

        for (int c = 0; c < GRAPH.getNumColors(); c++) {
            childColors.add(c);
        }

        childColors.removeAll(parentColors);
        childColors.removeAll(draftColors);
        childColors.remove(refColor);

        return childColors;
    }

    private boolean isSharedWithOnlyOneParent(CortexRecord cr, Set<Integer> parentColors, Set<Integer> draftColors) {
        int numParentsWithCoverage = 0;
        for (int pc : parentColors) {
            if (cr.getCoverage(pc) > 0) { numParentsWithCoverage++; }
        }

        int numDraftsWithCoverage = 0;
        for (int dc : draftColors) {
            if (cr.getCoverage(dc) > 0) { numDraftsWithCoverage++; }
        }

        return numParentsWithCoverage == 1 && numDraftsWithCoverage == 1;
    }

    private boolean isSharedWithSomeChildren(CortexRecord cr, Set<Integer> parentColors, Set<Integer> draftColors, int refColor) {
        Set<Integer> ignoreColors = new HashSet<>();
        ignoreColors.addAll(parentColors);
        ignoreColors.addAll(draftColors);
        ignoreColors.add(refColor);

        int numChildrenWithCoverage = 0;
        int numChildren = 0;
        for (int c = 0; c < cr.getNumColors(); c++) {
            if (!ignoreColors.contains(c)) {
                if (cr.getCoverage(c) > 0) {
                    numChildrenWithCoverage++;
                }
                numChildren++;
            }
        }

        return numChildrenWithCoverage > 1 && numChildrenWithCoverage < numChildren;
    }

    private boolean isSinglyConnected(CortexRecord cr) {
        boolean isSinglyConnected = true;
        for (int c = 0; c < cr.getNumColors(); c++) {
            if (cr.getCoverage(c) > 0) {
                if (!(cr.getInDegree(c) == 1 && cr.getOutDegree(c) == 1)) {
                    isSinglyConnected = false;
                }
            }
        }

        return isSinglyConnected;
    }

    private boolean hasUniqueCoordinates(CortexRecord cr, Set<Integer> draftColors) {
        int draftColor = -1;
        int numDraftsWithCoverage = 0;
        for (int dc : draftColors) {
            if (cr.getCoverage(dc) > 0) {
                draftColor = dc;
                numDraftsWithCoverage++;
            }
        }

        if (numDraftsWithCoverage == 1 && draftColor > -1) {
            KmerLookup kl = DRAFTS.get(GRAPH.getSampleName(draftColor));

            Set<Interval> its = kl.findKmer(cr.getKmerAsString());

            return its != null && its.size() == 1;
        }

        return false;
    }

    private boolean canWalkToCanonicalReference(CortexRecord cr, Set<Integer> childColors, int refColor) {
        TraversalEngine e = new TraversalEngineFactory()
                .graph(GRAPH)
                .make();

        for (int cc : childColors) {
            e.getConfiguration().setTraversalColor(cc);

            int found = 0;

            for (boolean goForward : Arrays.asList(true, false)) {
                e.setCursor(cr.getKmerAsString(), goForward);
                while (goForward ? e.hasNext() : e.hasPrevious()) {
                    CortexVertex cv = goForward ? e.next() : e.previous();

                    if (cv.getCr().getCoverage(refColor) > 0) {
                        found++;
                        break;
                    }
                }
            }

            if (found == 2) {
                return true;
            }
        }

        return false;
    }

    private List<Interval> getCanonicalReferenceCoordinates(CortexRecord cr, Set<Integer> childColors, int refColor) {
        TraversalEngine e = new TraversalEngineFactory()
                .graph(GRAPH)
                .make();

        for (int cc : childColors) {
            e.getConfiguration().setTraversalColor(cc);

            List<Interval> intervals = new ArrayList<>();

            for (boolean goForward : Arrays.asList(true, false)) {
                e.setCursor(cr.getKmerAsString(), goForward);
                while (goForward ? e.hasNext() : e.hasPrevious()) {
                    CortexVertex cv = goForward ? e.next() : e.previous();

                    if (cv.getCr().getCoverage(refColor) > 0) {
                        Set<Interval> its = REFERENCE.findKmer(cv.getSk());
                        if (its != null && its.size() == 1) {
                            intervals.add(its.iterator().next());
                        }

                        break;
                    }
                }
            }

            if (intervals.size() == 2) {
                Interval it0 = intervals.get(0);
                Interval it1 = intervals.get(1);

                if (it0.getContig().equals(it1.getContig())) {
                    return intervals;
                }
            }
        }

        return null;
    }

    private int getDraftColor(CortexRecord cr, Set<Integer> draftColors) {
        for (int dc : draftColors) {
            if (cr.getCoverage(dc) > 0) {
                return dc;
            }
        }

        return -1;
    }

    private Pair<String, String> callSimpleBubble(CortexRecord cr, int childColor, int bubbleColor) {
        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(childColor)
                .graph(GRAPH)
                .make();

        List<CortexVertex> childVertices = new ArrayList<>();
        childVertices.add(new CortexVertex(cr.getKmerAsString(), cr));

        for (boolean goForward : Arrays.asList(false, true)) {
            e.setCursor(cr.getKmerAsString(), goForward);
            while (goForward ? e.hasNext() : e.hasPrevious()) {
                CortexVertex cv = goForward ? e.next() : e.previous();

                if (goForward) { childVertices.add(cv); }
                else { childVertices.add(0, cv); }

                if (cv.getCr().getCoverage(bubbleColor) > 0) { break; }
            }
        }

        if (childVertices.get(0).getCr().getCoverage(bubbleColor) > 0) {
            e = new TraversalEngineFactory()
                    .traversalColor(bubbleColor)
                    .graph(GRAPH)
                    .make();

            List<CortexVertex> draftVertices = new ArrayList<>();
            draftVertices.add(childVertices.get(0));

            e.setCursor(childVertices.get(0).getSk(), true);
            while (e.hasNext()) {
                CortexVertex cv = e.next();

                draftVertices.add(cv);

                if (cv.getCr().getCoverage(childColor) > 0) { break; }
            }

            if (draftVertices.get(draftVertices.size() - 1).getSk().equals(childVertices.get(childVertices.size() - 1).getSk())) {
                return new Pair<>(TraversalEngine.toContig(childVertices), TraversalEngine.toContig(draftVertices));
            }
        }

        return null;
    }

    private Pair<String, String> trimToAlleles(Pair<String, String> haplotypes) {
        String s0 = haplotypes.getFirst();
        String s1 = haplotypes.getSecond();

        int s0start = 0, s0end = s0.length();
        int s1start = 0, s1end = s1.length();

        for (int i = 0, j = 0; i < s0.length() && j < s1.length(); i++, j++) {
            if (s0.charAt(i) != s1.charAt(j)) {
                s0start = i;
                s1start = j;
                break;
            }
        }

        for (int i = s0.length() - 1, j = s1.length() - 1; i >= 0 && j >= 0; i--, j--) {
            if (s0.charAt(i) != s1.charAt(j) || i == s0start - 1 || j == s1start - 1) {
                s0end = i + 1;
                s1end = j + 1;
                break;
            }
        }

        String[] pieces = new String[4];
        pieces[0] = s0.substring(0, s0start);
        pieces[1] = s0.substring(s0start, s0end);
        pieces[2] = s1.substring(s1start, s1end);
        pieces[3] = s0.substring(s0end, s0.length());

        return new Pair<>(pieces[1], pieces[2]);
    }

    private Interval getBubbleCanonicalCoordinates(String childHaplotype, int childColor, int refColor) {
        Interval is = null, ie = null;

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(childColor)
                .graph(GRAPH)
                .make();

        String firstKmer = childHaplotype.substring(0, GRAPH.getKmerSize());
        String lastKmer = childHaplotype.substring(childHaplotype.length() - GRAPH.getKmerSize(), childHaplotype.length());

        for (boolean goForward : Arrays.asList(false, true)) {
            e.setCursor(goForward ? lastKmer : firstKmer, goForward);
            while (goForward ? e.hasNext() : e.hasPrevious()) {
                CortexVertex cv = goForward ? e.next() : e.previous();

                if (cv.getCr().getCoverage(refColor) > 0) {
                    Set<Interval> intervals = REFERENCE.findKmer(cv.getSk());

                    if (intervals != null && intervals.size() == 1) {
                        if (goForward) {
                            ie = intervals.iterator().next();
                        } else {
                            is = intervals.iterator().next();
                        }

                        break;
                    }
                }
            }
        }

        if (is != null && ie != null) {
            return is.getStart() < ie.getStart() ? is : ie;
        }

        return null;
    }
}
