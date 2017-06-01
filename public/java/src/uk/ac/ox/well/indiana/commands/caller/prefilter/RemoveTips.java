package uk.ac.ox.well.indiana.commands.caller.prefilter;

import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultDirectedGraph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraphWriter;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.stoppingconditions.TipBeginningStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.TipEndStopper;
import uk.ac.ox.well.indiana.utils.traversal.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

@Description(text="Remove graph tips (chains of novel kmers only anchored at one end)")
public class RemoveTips extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Output
    public File out;

    @Output(fullName="tips_out", shortName="to", doc="Tips output file")
    public File tips_out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        Set<Integer> parentColors = new HashSet<>(GRAPH.getColorsForSampleNames(PARENTS));
        //Set<Integer> ignoreColors = new HashSet<>(GRAPH.getColorsForSampleNames(IGNORE));

        log.info("Colors:");
        log.info(" -   child: {}", GRAPH.getColorForSampleName(CHILD));
        log.info(" - parents: {}", GRAPH.getColorsForSampleNames(PARENTS));
        //log.info(" -  ignore: {}", GRAPH.getColorsForSampleNames(IGNORE));

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Finding tips")
                .message("records processed")
                .maxRecord(ROI.getNumRecords())
                .make(log);

        Set<CortexKmer> tips = new HashSet<>();
        int numTipChains = 0;

        for (CortexRecord rr : ROI) {
            if (!tips.contains(rr.getCortexKmer())) {
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfsToParents = null;
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfsToFree = null;

                for (boolean goForward : Arrays.asList(true, false)) {
                    dfsToParents = CortexUtils.dfs(GRAPH, rr.getKmerAsString(), childColor, parentColors, TipBeginningStopper.class, goForward);
                    dfsToFree = CortexUtils.dfs(GRAPH, rr.getKmerAsString(), childColor, parentColors, TipEndStopper.class, !goForward);

                    if (dfsToParents != null && dfsToParents.vertexSet().size() > 0 && dfsToFree != null && dfsToFree.vertexSet().size() > 0) {
                        break;
                    }
                }

                DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs = null;
                if (dfsToParents != null && dfsToParents.vertexSet().size() > 0 && dfsToFree != null && dfsToFree.vertexSet().size() > 0) {
                    dfs = new DefaultDirectedGraph<>(AnnotatedEdge.class);

                    Graphs.addGraph(dfs, dfsToParents);
                    Graphs.addGraph(dfs, dfsToFree);
                }


                if (dfs != null && dfs.vertexSet().size() > 0) {
                    numTipChains++;

                    log.debug("    tip chain {}, seed {}, {} vertices", numTipChains, rr.getKmerAsString(), dfs.vertexSet().size());

                    for (AnnotatedVertex av : dfs.vertexSet()) {
                        if (log.isDebugEnabled()) {
                            log.debug("    - {} {}", av.getKmer(), GRAPH.findRecord(new CortexKmer(av.getKmer())));
                        }

                        tips.add(new CortexKmer(av.getKmer()));
                    }
                }
            }

            pm.update();
        }

        log.info("Found {} tip kmer chains ({} kmers total)", numTipChains, tips.size());

        log.info("Writing...");

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(ROI.getHeader());

        CortexGraphWriter cgt = new CortexGraphWriter(tips_out);
        cgt.setHeader(ROI.getHeader());

        int numKept = 0, numExcluded = 0;
        for (CortexRecord rr : ROI) {
            if (!tips.contains(rr.getCortexKmer())) {
                cgw.addRecord(rr);
                numKept++;
            } else {
                cgt.addRecord(rr);
                numExcluded++;
            }
        }

        cgw.close();
        cgt.close();

        log.info("  {}/{} ({}%) kept, {}/{} ({}%) excluded",
                numKept,     ROI.getNumRecords(), 100.0f * (float) numKept / (float) ROI.getNumRecords(),
                numExcluded, ROI.getNumRecords(), 100.0f * (float) numExcluded / (float) ROI.getNumRecords()
        );
    }
}
