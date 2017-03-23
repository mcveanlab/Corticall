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
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedEdge;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedVertex;
import uk.ac.ox.well.indiana.utils.traversal.ChildTraversalStopper;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
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
        Set<CortexBinaryKmer> unused = new HashSet<>();

        int numFragments = 0;
        for (CortexRecord cr : NOVELS) {
            if (!seen.contains(cr.getCortexBinaryKmer())) {
                if (cr.getInDegree(0) <= 1 && cr.getOutDegree(0) <= 1) {
                    DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs = CortexUtils.dfs(GRAPHS, cr.getKmerAsString(), childColor, parentColors, ChildTraversalStopper.class);

                    int numNovel = 0;
                    //int totCov = 0;
                    for (AnnotatedVertex av : dfs.vertexSet()) {
                        CortexBinaryKmer cbk = new CortexBinaryKmer(av.getKmer().getBytes());
                        seen.add(cbk);

                        if (unused.contains(cbk)) { unused.remove(cbk); }

                        CortexRecord crn = GRAPHS.findRecord(new CortexKmer(av.getKmer()));
                        numNovel += CortexUtils.isNovelKmer(crn, childColor, parentColors) ? 1 : 0;

                        //totCov += crn.getCoverage(childColor);

                        /*
                        List<Integer> covs = new ArrayList<>();
                        List<String> edges = new ArrayList<>();

                        covs.add(crn.getCoverage(childColor));
                        edges.add(crn.getEdgesAsString(childColor));

                        for (int parentColor : parentColors) {
                            covs.add(crn.getCoverage(parentColor));
                            edges.add(crn.getEdgesAsString(parentColor));
                        }

                        log.info("    {} {} {}", crn.getCortexKmer(), Joiner.on(' ').join(covs), Joiner.on(' ').join(edges));
                        */
                    }

                    log.info("    fragment {}: {} {} {}", numFragments, dfs.vertexSet().size(), dfs.edgeSet().size(), numNovel);

                    numFragments++;
                } else {
                    unused.add(cr.getCortexBinaryKmer());
                }
            }

            pm.update("records processed (" + numFragments + " fragments constructed so far, " + unused.size() + " unused)");
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
