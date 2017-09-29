package uk.ac.ox.well.cortexjdk.playground.assemblies;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by kiran on 26/08/2017.
 */
public class AssembleWithReferenceAssist extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="links", shortName="l", doc="Links")
    public CortexLinks LINKS;

    @Argument(fullName="ref", shortName="R", doc="Reference")
    public KmerLookup REF;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TraversalEngine e0 = new TraversalEngineFactory()
                .traversalColor(0)
                .graph(GRAPH)
                .make();

        TraversalEngine el = new TraversalEngineFactory()
                .traversalColor(0)
                .graph(GRAPH)
                .links(LINKS)
                //.references(REF)
                .make();

        for (CortexRecord cr : GRAPH) {
            if (cr.getCoverage(0) > 0) {
                List<CortexVertex> w0 = e0.walk(cr.getKmerAsString());

                if (w0.size() > 100) {
                    List<CortexVertex> wl = el.walk(cr.getKmerAsString());

                    if (wl.size() >= w0.size() + 200) {
                        log.info("{} {} {}", w0.size(), wl.size(), cr);
                        log.info("{}", wl.get(0));
                        log.info("{}", wl.get(wl.size() - 1));

                        for (CortexVertex cv : wl) {
                            log.info("{} {} {}", cv.getSk(), cv.getSources(), REF.find(cv.getSk()));
                        }
                    }
                }
            }
        }
    }

    private List<CortexVertex> annotateContig(List<CortexVertex> cvs) {
        List<CortexVertex> cvas = new ArrayList<>();

        return null;
    }
}
