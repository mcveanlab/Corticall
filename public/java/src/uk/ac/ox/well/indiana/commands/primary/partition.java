package uk.ac.ox.well.indiana.commands.primary;

import com.google.api.client.util.Joiner;
import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexBinaryKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.traversal.AbstractTraversalStopper;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedEdge;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedVertex;
import uk.ac.ox.well.indiana.utils.traversal.ChildTraversalStopper;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class partition extends Module {
    @Argument(fullName="novels", shortName="n", doc="Novels")
    public CortexGraph NOVELS;

    @Output
    public File out;

    @Override
    public void execute() {
        CortexGraph cg = new CortexGraph(NOVELS.getCortexFile());

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing graph...")
                .message("records processed")
                .updateRecord(NOVELS.getNumRecords() / 10)
                .maxRecord(NOVELS.getNumRecords())
                .make(log);

        Set<CortexBinaryKmer> seen = new HashSet<>();
        Set<CortexBinaryKmer> unused = new HashSet<>();

        int numFragments = 0;
        for (CortexRecord cr : NOVELS) {
            if (!seen.contains(cr.getCortexBinaryKmer())) {
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs = CortexUtils.dfs(cg, cr.getKmerAsString(), 0, null, NovelKmerAggregator.class);

                for (AnnotatedVertex av : dfs.vertexSet()) {
                    CortexBinaryKmer cbk = new CortexBinaryKmer(av.getKmer().getBytes());

                    seen.add(cbk);
                }

                numFragments++;
            }

            pm.update("records processed (" + numFragments + " fragments constructed so far, " + unused.size() + " unused)");
        }
    }

    class NovelKmerAggregator extends AbstractTraversalStopper<AnnotatedVertex, AnnotatedEdge> {
        @Override
        public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions, int size, int edges, Set<Integer> childColors, Set<Integer> parentColors) {
            return cr.getInDegree(0) == 0 || cr.getOutDegree(0) == 0;
        }

        @Override
        public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions, int size, int edges, Set<Integer> childColors, Set<Integer> parentColors) {
            return false;
        }

        @Override
        public int maxJunctionsAllowed() {
            return 0;
        }
    }
}
