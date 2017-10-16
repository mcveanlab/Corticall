package uk.ac.ox.well.cortexjdk.commands.call.call;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.ConnectivityAnnotations;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.BubbleOpeningStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.AND;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

/**
 * Created by kiran on 30/08/2017.
 */
public class CallBreakpoints extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="links", shortName="l", doc="Links")
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName="references", shortName="R", doc="References")
    public HashMap<String, KmerLookup> REFERENCES;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CortexKmer, Boolean> seen = new HashMap<>();
        /*
        for (CortexRecord cr : ROI) {
            seen.put(cr.getCortexKmer(), false);
        }
        */
        seen.put(new CortexKmer("AACCCAAAAAAGTATATCCTTAACCCTGAGAGATGGGAACCCTAAAC"), false);

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(GRAPH.getColorForSampleName(ROI.getSampleName(0)))
                .joiningColors(GRAPH.getColorsForSampleNames(REFERENCES.keySet()))
                .graph(GRAPH)
                .links(LINKS)
                .references(REFERENCES.values())
                .make();

        for (CortexKmer ck : seen.keySet()) {
            if (!seen.get(ck)) {
                List<CortexVertex> w = e.walk(ck.getKmerAsString());

                log.info("{}", ck);
                for (CortexVertex v : w) {
                    log.info("\t{} {}", v, REFERENCES.get("PG0051-C.ERR019061").find(v.getSk()));

                    if (seen.containsKey(v.getCk())) {
                        seen.put(v.getCk(), true);
                    }
                }
            }
        }
    }
}
