package uk.ac.ox.well.indiana.commands.primary;

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
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedEdge;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedVertex;
import uk.ac.ox.well.indiana.utils.traversal.ChildTraversalStopper;

import java.io.File;
import java.util.HashSet;
import java.util.Set;

public class partition extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPHS;

    @Argument(fullName="novels", shortName="n", doc="Novels")
    public CortexGraph NOVELS;

    @Argument(fullName="parent", shortName="p", doc="Parent")
    public HashSet<String> PARENTS;

    @Output
    public File out;

    @Override
    public void execute() {
        int childColor = GRAPHS.getColorForSampleName(NOVELS.getSampleName(0));
        Set<Integer> parentColors = getParentColors();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing graph...")
                .message("records processed")
                .updateRecord(NOVELS.getNumRecords() / 10)
                .maxRecord(NOVELS.getNumRecords())
                .make(log);

        Set<CortexBinaryKmer> seen = new HashSet<>();

        int numFragments = 0;
        for (CortexRecord cr : NOVELS) {
            if (!seen.contains(cr.getCortexBinaryKmer())) {
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs = CortexUtils.dfs(GRAPHS, cr.getKmerAsString(), childColor, parentColors, ChildTraversalStopper.class);

                log.info("  fragment");

                int numNovel = 0;
                for (AnnotatedVertex av : dfs.vertexSet()) {
                    seen.add(new CortexBinaryKmer(av.getKmer().getBytes()));

                    CortexRecord crn = GRAPHS.findRecord(new CortexKmer(av.getKmer()));
                    numNovel += CortexUtils.isNovelKmer(crn, childColor, parentColors) ? 1 : 0;

                    log.info("    {}", crn);
                }

                log.info("    fragment {}: {} {} {}", numFragments, dfs.vertexSet().size(), dfs.edgeSet().size(), numNovel);

                numFragments++;
            }

            pm.update("records processed (" + numFragments + " constructed so far)");
        }
    }

    private Set<Integer> getParentColors() {
        Set<Integer> parentColors = new HashSet<>();

        for (String parent : PARENTS) {
            parentColors.add(GRAPHS.getColorForSampleName(parent));
        }

        return parentColors;
    }
}
