package uk.ac.ox.well.cortexjdk.commands.prefilter;

import org.jetbrains.annotations.NotNull;
import org.jgrapht.Graph;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Description;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.*;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContaminantStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

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

    @Output(fullName="removed_out", shortName="ro", doc="Contam output file")
    public File contam_out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        Set<Integer> parentColors = new HashSet<>(GRAPH.getColorsForSampleNames(PARENTS));

        log.info("Colors:");
        log.info(" -   child: {}", GRAPH.getColorForSampleName(CHILD));
        log.info(" - parents: {}", GRAPH.getColorsForSampleNames(PARENTS));

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Finding contamination")
                .message("records processed")
                .maxRecord(CONTAM.getNumRecords())
                .make(log);

        Set<CanonicalKmer> roiKmers = new HashSet<>();
        for (CortexRecord rc : ROI) {
            roiKmers.add(rc.getCanonicalKmer());
        }

        Set<CanonicalKmer> contamKmers = new HashSet<>();
        int numContamChains = 0;

        TraversalEngine e = new TraversalEngineFactory()
                .traversalDirection(BOTH)
                .combinationOperator(OR)
                .traversalColor(childColor)
                .joiningColors(parentColors)
                .stoppingRule(ContaminantStopper.class)
                .rois(ROI)
                .graph(GRAPH)
                .make();

        for (CortexRecord cr : CONTAM) {
            if (roiKmers.contains(cr.getCanonicalKmer()) && !contamKmers.contains(cr.getCanonicalKmer())) {
                Graph<CortexVertex, CortexEdge> dfs = e.dfs(cr.getKmerAsString());

                if (dfs.vertexSet().size() > 0) {
                    for (CortexVertex av : dfs.vertexSet()) {
                        CanonicalKmer ck = av.getCanonicalKmer();

                        contamKmers.add(ck);
                    }

                    numContamChains++;
                }
            }

            pm.update();
        }

        log.info("Found {} contamination kmer chains ({} kmers total)", numContamChains, contamKmers.size());

        log.info("Writing...");

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(makeCortexHeader());

        CortexGraphWriter cgc = new CortexGraphWriter(contam_out);
        cgc.setHeader(ROI.getHeader());

        int numKept = 0, numExcluded = 0;
        for (CortexRecord rr : ROI) {
            if (!contamKmers.contains(rr.getCanonicalKmer())) {
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
