package uk.ac.ox.well.indiana.commands.inheritance;

import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.stoppingconditions.BubbleOpeningStopper;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngineFactory;

import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.AND;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

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

        for (CortexRecord cr : GRAPH) {
            if (isSharedWithOnlyOneParent(cr, parentColors, draftColors) &&
                isSharedWithSomeChildren(cr, parentColors, draftColors, refColor) &&
                isSinglyConnected(cr) &&
                hasUniqueCoordinates(cr, draftColors) &&
                canWalkToCanonicalReference(cr, childColors, refColor)) {

                int draftColor = getDraftColor(cr, draftColors);

                List<Interval> intervals = getCanonicalReferenceCoordinates(cr, childColors, refColor);

                if (intervals != null) {
                    log.info("{} {} {} {}", draftColor, GRAPH.getSampleName(draftColor), intervals, cr);
                    for (String id : DRAFTS.keySet()) {
                        log.info("  {} {}", id, DRAFTS.get(id).findKmer(cr.getKmerAsString()));
                    }
                }
            }
        }
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
}
