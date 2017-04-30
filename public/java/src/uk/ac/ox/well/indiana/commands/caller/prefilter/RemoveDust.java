package uk.ac.ox.well.indiana.commands.caller.prefilter;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraphWriter;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.stoppingconditions.DustStopper;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedEdge;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedVertex;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

public class RemoveDust extends Module {
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

    @Output(fullName="dust_out", shortName="do", doc="Dust output file")
    public File dust_out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        Set<Integer> parentColors = new HashSet<>(GRAPH.getColorsForSampleNames(PARENTS));

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Finding dust")
                .message("records processed")
                .maxRecord(ROI.getNumRecords())
                .make(log);

        Set<CortexKmer> dust = new HashSet<>();
        int numDustChains = 0;

        for (CortexRecord rr : ROI) {
            if (!dust.contains(rr.getCortexKmer()) && CortexUtils.isLowComplexity(rr, 0)) {
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs = CortexUtils.dfs_and(GRAPH, rr.getKmerAsString(), childColor, parentColors, DustStopper.class);

                if (dfs != null && dfs.vertexSet().size() > 0) {
                    numDustChains++;

                    log.debug("    dust chain {}, seed {}, {} vertices", numDustChains, rr.getKmerAsString(), dfs.vertexSet().size());

                    for (AnnotatedVertex av : dfs.vertexSet()) {
                        if (log.isDebugEnabled()) {
                            log.debug("    - {} {}", av.getKmer(), GRAPH.findRecord(new CortexKmer(av.getKmer())));
                        }

                        dust.add(new CortexKmer(av.getKmer()));
                    }

                }
            }

            pm.update();
        }

        log.info("Found {} dust kmer chains ({} kmers total)", numDustChains, dust.size());

        log.info("Writing...");

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(ROI.getHeader());

        CortexGraphWriter cgo = new CortexGraphWriter(dust_out);
        cgo.setHeader(ROI.getHeader());

        int numKept = 0, numExcluded = 0;
        for (CortexRecord rr : ROI) {
            if (!dust.contains(rr.getCortexKmer())) {
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
                numKept,     ROI.getNumRecords(), numKept / ROI.getNumRecords(),
                numExcluded, ROI.getNumRecords(), numExcluded / ROI.getNumRecords()
        );
    }
}
