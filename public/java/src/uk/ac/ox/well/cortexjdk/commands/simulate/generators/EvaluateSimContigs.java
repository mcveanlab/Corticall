package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import com.google.common.base.Joiner;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.util.*;

public class EvaluateSimContigs extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="links", shortName="l", doc="Links")
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(ROI.getColor(0).getSampleName());

        TraversalEngine el = new TraversalEngineFactory()
                .traversalColor(childColor)
                .graph(GRAPH)
                .links(LINKS)
                .make();

        TraversalEngine er = new TraversalEngineFactory()
                .traversalColor(childColor)
                .graph(GRAPH)
                .make();

        Map<CanonicalKmer, String> seen = new HashMap<>();
        for (CortexRecord rr : ROI) {
            seen.put(rr.getCanonicalKmer(), null);
        }

        for (CortexRecord rr : ROI) {
            if (seen.containsKey(rr.getCanonicalKmer())) {
                log.info("{}", seen.get(rr.getCanonicalKmer()));
            } else {
                List<CortexVertex> erw = er.walk(rr.getKmerAsString());
                List<CortexVertex> elw = el.walk(rr.getKmerAsString());

                String out = Joiner.on(" ").join(rr.getKmerAsString(), erw.size(), elw.size());

                log.info("{}", out);

                seen.put(rr.getCanonicalKmer(), out);
                for (CortexVertex cv : erw) {
                    if (seen.containsKey(cv.getCanonicalKmer()) && seen.get(cv.getCanonicalKmer()) == null) {
                        seen.put(cv.getCanonicalKmer(), out);
                    }
                }
            }
        }
    }
}
