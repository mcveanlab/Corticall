package uk.ac.ox.well.indiana.commands.guided;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksMap;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngineFactory;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by kiran on 03/08/2017.
 */
public class EmitGuidedAssembly extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="links", shortName="l", doc="Links")
    public HashMap<CortexLinks, String> LINKS;

    @Argument(fullName="refs", shortName="R", doc="References")
    public HashMap<KmerLookup, String> REFERENCES;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="roi", shortName="r", doc="ROIs")
    public CortexGraph ROI;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(new ArrayList<>(REFERENCES.values()));

        log.info("Colors:");
        log.info("  -   child: {} {}", CHILD, childColor);
        log.info("  - parents: {} {}", REFERENCES.values(), parentColors);

        TraversalEngine e = new TraversalEngineFactory()
                .graph(GRAPH)
                .links(LINKS)
                .traversalColor(childColor)
                .rois(ROI)
                .make();

        log.info("Processing contigs:");

        for (KmerLookup REF : REFERENCES.keySet()) {
            ReferenceSequence rseq;
            while ((rseq = REF.getReferenceSequence().nextSequence()) != null) {
                log.info("  {}", rseq.getName());

                String seq = rseq.getBaseString();

                Map<String, Boolean> signalKmers = new HashMap<>();

                for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                    String sk = seq.substring(i, i + GRAPH.getKmerSize());
                    CortexKmer ck = new CortexKmer(sk);
                    CortexRecord cr = GRAPH.findRecord(ck);
                    Set<Interval> its = REF.findKmer(sk);

                    if (isSignalKmer(childColor, parentColors, cr, its)) {
                        signalKmers.put(sk, false);
                    }
                }

                log.info("  - {} {}", seq.length() - GRAPH.getKmerSize(), signalKmers.size());

                for (String signalKmer : signalKmers.keySet()) {
                    List<CortexVertex> cvs = new ArrayList<>();

                    for (boolean goForward : Arrays.asList(true, false)) {
                        e.setCursor(signalKmer, goForward);
                        while (goForward ? e.hasNext() : e.hasPrevious()) {
                            CortexVertex cv = goForward ? e.next() : e.previous();

                            if (goForward) {
                                cvs.add(cv);
                            } else {
                                cvs.add(0, cv);
                            }

                            log.info("{} {} {} {} {}", signalKmer, goForward, cv.getSk(), cv.getSources(), REF.findKmer(cv.getSk()));
                        }
                    }
                }
            }
        }
    }

    private boolean isSignalKmer(int childColor, List<Integer> parentColors, CortexRecord cr, Set<Interval> its) {
        return cr != null &&
               cr.getCoverage(childColor) > 0 &&
               cr.getInDegree(childColor) == 1 &&
               cr.getOutDegree(childColor) == 1 &&
               hasNoLinks(cr) &&
               numParentsWithCoverage(cr, parentColors) == 1 &&
               singlyConnected(cr, parentColors) &&
               its.size() == 1;
    }

    private boolean hasNoLinks(CortexRecord cr) {
        for (CortexLinks lm : LINKS.keySet()) {
            if (lm.containsKey(cr.getCortexKmer())) {
                return false;
            }
        }

        return true;
    }

    private int numParentsWithCoverage(CortexRecord cr, List<Integer> parentColors) {
        int numParentsWithCoverage = 0;

        for (int c : parentColors) {
            if (cr.getCoverage(c) > 0) {
                numParentsWithCoverage++;
            }
        }

        return numParentsWithCoverage;
    }

    private boolean singlyConnected(CortexRecord cr, int c) {
        return cr.getCoverage(c) > 0 && cr.getInDegree(c) == 1 && cr.getOutDegree(c) == 1;
    }

    private boolean singlyConnected(CortexRecord cr, List<Integer> parentColors) {
        for (int c : parentColors) {
            if (singlyConnected(cr, c)) {
                return true;
            }
        }

        return false;
    }
}
