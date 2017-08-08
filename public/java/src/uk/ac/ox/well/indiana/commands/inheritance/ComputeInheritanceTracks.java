package uk.ac.ox.well.indiana.commands.inheritance;

import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.io.PrintStream;
import java.util.*;

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

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Set<Integer> parentColors = new TreeSet<>(GRAPH.getColorsForSampleNames(PARENTS));
        Set<Integer> draftColors = new TreeSet<>(GRAPH.getColorsForSampleNames(new ArrayList<>(DRAFTS.keySet())));
        int refColor = GRAPH.getColorForSampleName("ref");

        for (CortexRecord cr : GRAPH) {
            if (isSharedWithOnlyOneParent(cr, parentColors, draftColors) &&
                isSharedWithSomeChildren(cr, parentColors, draftColors, refColor) &&
                isSinglyConnected(cr) &&
                hasUniqueCoordinates(cr, draftColors)) {
                log.info("{}", cr);
                for (String id : DRAFTS.keySet()) {
                    log.info("  {} {}", id, DRAFTS.get(id).findKmer(cr.getKmerAsString()));
                }
            }
        }
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
}
