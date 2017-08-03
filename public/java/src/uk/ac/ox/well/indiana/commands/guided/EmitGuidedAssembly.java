package uk.ac.ox.well.indiana.commands.guided;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
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

    //@Argument(fullName="links", shortName="l", doc="Links")
    //public HashMap<CortexLinks, String> LINKS;

    @Argument(fullName="refs", shortName="R", doc="References")
    public HashMap<IndexedFastaSequenceFile, String> REFERENCES;

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
                //.links(LINKS)
                .traversalColor(childColor)
                .rois(ROI)
                .make();

        for (IndexedFastaSequenceFile REF : REFERENCES.keySet()) {
            ReferenceSequence rseq;
            while ((rseq = REF.nextSequence()) != null) {
                String seq = rseq.getBaseString();

                Map<String, Boolean> signalKmers = new HashMap<>();

                for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                    String sk = seq.substring(i, i + GRAPH.getKmerSize());
                    CortexKmer ck = new CortexKmer(sk);
                    CortexRecord cr = GRAPH.findRecord(ck);

                    if (cr != null &&
                        cr.getCoverage(childColor) > 0 &&
                        cr.getInDegree(childColor) == 1 && cr.getOutDegree(childColor) == 1 &&
                        numParentsWithCoverage(cr, parentColors) == 1 &&
                        singlyConnected(cr, parentColors)) {

                        signalKmers.put(sk, false);
                    }
                }

                log.info("{} {} {}", rseq.getName(), seq.length() - GRAPH.getKmerSize(), signalKmers.size());
            }
        }
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
