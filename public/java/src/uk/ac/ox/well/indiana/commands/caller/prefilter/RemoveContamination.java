package uk.ac.ox.well.indiana.commands.caller.prefilter;

import org.jetbrains.annotations.NotNull;
import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.*;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedEdge;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedVertex;
import uk.ac.ox.well.indiana.utils.stoppingconditions.ContaminantStopper;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

@Description(text="Remove chains of contaminating kmers")
public class RemoveContamination extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="contamination", shortName="contam", doc="Contam")
    public CortexGraph CONTAM;

    @Output
    public File out;

    @Output(fullName="contam_out", shortName="co", doc="Contam output file")
    public File contam_out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        Set<Integer> parentColors = new HashSet<>(GRAPH.getColorsForSampleNames(PARENTS));

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing records")
                .message("records processed")
                .maxRecord(CONTAM.getNumRecords())
                .make(log);

        Set<CortexKmer> roiKmers = new HashSet<>();
        for (CortexRecord rc : ROI) {
            roiKmers.add(rc.getCortexKmer());
        }

        Set<CortexKmer> contamKmers = new HashSet<>();
        int numContamChains = 0;

        for (CortexRecord cr : CONTAM) {
            pm.update();

            if (roiKmers.contains(cr.getCortexKmer()) && !contamKmers.contains(cr.getCortexKmer())) {
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs = CortexUtils.dfs(GRAPH, cr.getKmerAsString(), childColor, parentColors, ContaminantStopper.class);

                if (dfs.vertexSet().size() > 0) {
                    for (AnnotatedVertex av : dfs.vertexSet()) {
                        CortexKmer ck = new CortexKmer(av.getKmer());

                        contamKmers.add(ck);
                    }

                    numContamChains++;
                }
            }
        }

        log.info("Found {} contamination kmer chains ({} kmers total)", numContamChains, contamKmers.size());

        log.info("Writing...");

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(makeCortexHeader());

        CortexGraphWriter cgc = new CortexGraphWriter(contam_out);
        cgc.setHeader(ROI.getHeader());

        int numKept = 0, numExcluded = 0;
        for (CortexRecord rr : ROI) {
            if (!contamKmers.contains(rr.getCortexKmer())) {
                cgw.addRecord(rr);
                numKept++;
            } else {
                cgc.addRecord(rr);
                numExcluded++;
            }
        }

        cgw.close();
        cgc.close();

        log.info("  {}/{} ({}%) kept, {}/{} ({}%) excluded",
                numKept,     ROI.getNumRecords(), 100.0f * (float) numKept / (float) ROI.getNumRecords(),
                numExcluded, ROI.getNumRecords(), 100.0f * (float) numExcluded / (float) ROI.getNumRecords()
        );
    }

    @NotNull
    private CortexHeader makeCortexHeader() {
        CortexHeader ch = new CortexHeader();
        ch.setVersion(6);
        ch.setNumColors(1);
        ch.setKmerSize(GRAPH.getKmerSize());
        ch.setKmerBits(GRAPH.getKmerBits());

        CortexColor cc = new CortexColor();
        cc.setCleanedAgainstGraph(false);
        cc.setCleanedAgainstGraphName("");
        cc.setErrorRate(0);
        cc.setLowCovgKmersRemoved(false);
        cc.setLowCovgSupernodesRemoved(false);
        cc.setTipClippingApplied(false);
        cc.setTotalSequence(0);
        cc.setSampleName(CHILD);

        ch.addColor(cc);

        return ch;
    }
}
