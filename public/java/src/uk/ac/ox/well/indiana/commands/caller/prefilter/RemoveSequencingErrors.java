package uk.ac.ox.well.indiana.commands.caller.prefilter;

import com.google.common.collect.Lists;
import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.stoppingconditions.ChildTraversalStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.UniquePathStopper;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedEdge;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngineFactory;

import java.util.*;

public class RemoveSequencingErrors extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="roi", shortName="r", doc="ROIs")
    public CortexGraph ROI;

    @Override
    public void execute() {
        TraversalEngine ce = new TraversalEngineFactory()
                .graph(GRAPH)
                .traversalSamples(CHILD)
                .joiningSamples(PARENTS)
                .stopper(new UniquePathStopper())
                .make();

        TraversalEngine pe = new TraversalEngineFactory()
                .graph(GRAPH)
                .traversalSamples(CHILD)
                .joiningSamples(PARENTS)
                //.stopper(new SharedPathStopper())
                .make();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Finding shared kmers")
                .message("records processed")
                .maxRecord(ROI.getNumRecords())
                .make(log);

        Set<CortexKmer> errorKmers = new HashSet<>();

        for (CortexRecord rr : ROI) {
            if (!errorKmers.contains(rr.getCortexKmer())) {
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfsUnique = ce.dfs(rr.getKmerAsString());
                //DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfsShared = pe.dfs(dfsUnique);



                //DirectedGraph<AnnotatedVertex, AnnotatedEdge> gchild = CortexUtils.dfs(GRAPH, rr.getKmerAsString(), childColor, parentColors, ChildTraversalStopper.class);
                //DirectedGraph<AnnotatedVertex, AnnotatedEdge> gparent1 = CortexUtils.dfs(GRAPH, null, rr.getKmerAsString(), )
            }

            /*
            if (!sharedKmers.contains(rr.getCortexKmer())) {
                CortexRecord cr = GRAPH.findRecord(rr.getCortexKmer());

                for (int c = 0; c < GRAPH.getNumColors(); c++) {
                    if (c != childColor && !parentColors.contains(c) && !ignoreColors.contains(c) && cr.getCoverage(c) > 0) {
                        sharedKmers.add(rr.getCortexKmer());

                        log.debug("{}", cr);

                        break;
                    }
                }
            }
            */

            pm.update();
        }

        log.info("Found {} shared kmers", errorKmers.size());

        log.info("Writing...");

        /*
        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(ROI.getHeader());

        CortexGraphWriter cgo = new CortexGraphWriter(shared_out);
        cgo.setHeader(ROI.getHeader());

        int numKept = 0, numExcluded = 0;
        for (CortexRecord rr : ROI) {
            if (!sharedKmers.contains(rr.getCortexKmer())) {
                cgw.addRecord(rr);
                numKept++;
            } else {
                cgo.addRecord(rr);
                numExcluded++;
            }
        }

        cgw.close();
        cgo.close();

        log.info("  {}/{} ({}%) kept, {}/{} ({}%) excluded",
                numKept,     ROI.getNumRecords(), 100.0f * (float) numKept / (float) ROI.getNumRecords(),
                numExcluded, ROI.getNumRecords(), 100.0f * (float) numExcluded / (float) ROI.getNumRecords()
        );
        */
    }
}
