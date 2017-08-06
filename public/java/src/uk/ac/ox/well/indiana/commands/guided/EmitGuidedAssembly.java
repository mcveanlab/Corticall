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
import java.sql.Ref;
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

        Map<String, KmerLookup> ids = new TreeMap<>();
        for (KmerLookup kl : REFERENCES.keySet()) {
            ids.put(REFERENCES.get(kl), kl);
        }

        for (String id : ids.keySet()) {
            KmerLookup ref = ids.get(id);

            ReferenceSequence rseq;
            while ((rseq = ref.getReferenceSequence().nextSequence()) != null) {
                log.info("  {}", rseq.getName());

                String seq = rseq.getBaseString();

                Map<String, Boolean> signalKmers = new HashMap<>();

                for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                    String sk = seq.substring(i, i + GRAPH.getKmerSize());
                    CortexKmer ck = new CortexKmer(sk);
                    CortexRecord cr = GRAPH.findRecord(ck);
                    Set<Interval> its = ref.findKmer(sk);

                    if (isSignalKmer(childColor, parentColors, cr, its, rseq)) {
                        signalKmers.put(sk, false);
                    }
                }

                log.info("  - {} {}", seq.length() - GRAPH.getKmerSize(), signalKmers.size());

                for (String signalKmer : signalKmers.keySet()) {
                    List<CortexVertex> cvs = new ArrayList<>();

                    CortexVertex cv0 = new CortexVertex(signalKmer,
                                                        GRAPH.findRecord(new CortexKmer(signalKmer)),
                                                        ref.findKmer(signalKmer).iterator().next(),
                                                        Collections.singleton(REFERENCES.get(ref)));
                    cvs.add(cv0);

                    for (boolean goForward : Arrays.asList(true, false)) {
                        Interval lastConfidentInterval = cv0.getLocus();
                        boolean isAscending = true;

                        e.setCursor(signalKmer, goForward);
                        while (goForward ? e.hasNext() : e.hasPrevious()) {
                            CortexVertex cv = goForward ? e.next() : e.previous();

                            Set<Interval> its = ref.findKmer(cv.getSk());
                            Interval it = selectInterval(lastConfidentInterval, its);

                            if (goForward) {
                                cvs.add(cv);
                            } else {
                                cvs.add(0, cv);
                            }

                            log.info("{} {} {} {} {}", goForward, cv.getSk(), cv.getSources(), it, ref.findKmer(cv.getSk()));

                            if (goForward ? !e.hasNext() : !e.hasPrevious()) {
                                log.info("  {}", cv.getCr());
                            }
                        }
                    }
                }
            }
        }
    }

    private Interval selectInterval(Interval lastConfidentInterval, Set<Interval> its) {
        if (lastConfidentInterval != null) {
            for (Interval it : its) {
                if (it.intersects(lastConfidentInterval)) {
                    return it;
                }
            }

            return null;
        }

        return its.iterator().next();
    }

    private boolean isSignalKmer(int childColor, List<Integer> parentColors, CortexRecord cr, Set<Interval> its, ReferenceSequence rseq) {
        return cr != null &&
               cr.getCoverage(childColor) > 0 &&
               cr.getInDegree(childColor) == 1 &&
               cr.getOutDegree(childColor) == 1 &&
               hasNoLinks(cr) &&
               numParentsWithCoverage(cr, parentColors) == 1 &&
               singlyConnected(cr, parentColors) &&
               its.size() == 1 &&
               rseq.getName().contains(its.iterator().next().getContig());
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
