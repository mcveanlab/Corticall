package uk.ac.ox.well.indiana.commands.playground.caller;

import org.jetbrains.annotations.NotNull;
import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.*;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedEdge;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedVertex;
import uk.ac.ox.well.indiana.utils.traversal.ContaminantStopper;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

public class RemoveContamination extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="contamination", shortName="contam", doc="Contam")
    public CortexGraph CONTAM;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Output
    public File out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        Set<Integer> parentColors = new HashSet<>();
        for (String parent : PARENTS) {
            parentColors.add(GRAPH.getColorForSampleName(parent));
        }

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing records")
                .message("records processed")
                .maxRecord(CONTAM.getNumRecords())
                .make(log);

        Set<CortexKmer> seen = new HashSet<>();
        for (CortexRecord contam : CONTAM) {
            pm.update();

            if (!seen.contains(contam.getCortexKmer())) {
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs = CortexUtils.dfs(GRAPH, contam.getKmerAsString(), childColor, parentColors, ContaminantStopper.class);

                for (AnnotatedVertex av : dfs.vertexSet()) {
                    CortexKmer ck = new CortexKmer(av.getKmer());

                    seen.add(ck);
                }
            }
        }

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(makeCortexHeader());

        for (CortexRecord cr : ROI) {
            if (!seen.contains(cr.getCortexKmer())) {
                cgw.addRecord(cr);
            }
        }

        cgw.close();
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
