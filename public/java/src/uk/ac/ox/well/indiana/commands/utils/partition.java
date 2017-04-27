package uk.ac.ox.well.indiana.commands.utils;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexBinaryKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.traversal.*;

import java.io.File;
import java.util.*;

public class partition extends Module {
    @Argument(fullName="graphs", shortName="g", doc="Graphs")
    public CortexGraph GRAPHS;

    @Argument(fullName="novels", shortName="n", doc="Novels")
    public CortexGraph NOVELS;

    @Output
    public File out;

    @Override
    public void execute() {
        /*
        DB db = DBMaker
                .fileDB(out)
                .make();
                */

        int childColor = GRAPHS.getColorForSampleName(NOVELS.getSampleName(0));

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing graph...")
                .message("records processed")
                .updateRecord(NOVELS.getNumRecords() / 10)
                .maxRecord(NOVELS.getNumRecords())
                .make(log);

        Set<CortexBinaryKmer> seen = new HashSet<>();

        /*
        NavigableSet partition = db.treeSet("partition")
                .counterEnable()
                .createOrOpen();
                */

        Set<DirectedGraph<AnnotatedVertex, AnnotatedEdge>> partition = new HashSet<>();

        for (CortexRecord cr : NOVELS) {
            pm.update("records processed (" + partition.size() + " fragments constructed so far)");

            if (!seen.contains(cr.getCortexBinaryKmer())) {
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs = CortexUtils.dfs(GRAPHS, cr.getKmerAsString(), childColor, null, NovelKmerAggregationStopper.class);

                for (AnnotatedVertex av : dfs.vertexSet()) {
                    CortexBinaryKmer cbk = new CortexBinaryKmer(av.getKmer().getBytes());

                    seen.add(cbk);
                }

                partition.add(dfs);

                log.info("  {} {}", dfs.vertexSet().size(), dfs.edgeSet().size());
                log.info("");
            }
        }

        log.info("  Fragments: {}", partition.size());

        /*
        db.commit();

        db.close();
        */
    }
}
