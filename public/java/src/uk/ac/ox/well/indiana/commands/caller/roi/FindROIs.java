package uk.ac.ox.well.indiana.commands.caller.roi;

import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.*;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

@Description(text="Generate a (very liberal) list of kmers that identify potential de novo mutations")
public class FindROIs extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Output
    public File out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = new ArrayList<>(GRAPH.getColorsForSampleNames(PARENTS));

        log.info("Colors:");
        log.info(" -   child: {}", GRAPH.getColorForSampleName(CHILD));
        log.info(" - parents: {}", GRAPH.getColorsForSampleNames(PARENTS));

        log.info("Color: {} {}", CHILD, childColor);

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing graph...")
                .message("records processed")
                .maxRecord(GRAPH.getNumRecords())
                .make(log);

        long numNovelRecords = 0L;

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(makeCortexHeader(childColor));

        for (CortexRecord cr : GRAPH) {
            if (isNovel(cr, parentColors, childColor)) {
                CortexRecord novelCr = new CortexRecord(
                    cr.getBinaryKmer(),
                    new int[] { cr.getCoverages()[childColor] },
                    new byte[] { cr.getEdges()[childColor] },
                    cr.getKmerSize(), cr.getKmerBits()
                );

                cgw.addRecord(novelCr);

                numNovelRecords++;
            }

            pm.update("records processed (" + numNovelRecords + " novel so far)");
        }

        cgw.close();
    }

    private boolean isNovel(CortexRecord cr, List<Integer> parentColors, int childColor) {
        boolean parentsLackCoverage = true;

        for (int c : parentColors) {
            parentsLackCoverage &= cr.getCoverage(c) == 0;
        }

        boolean childHasCoverage = cr.getCoverage(childColor) > 0;

        return childHasCoverage && parentsLackCoverage;
    }

    @NotNull
    private CortexHeader makeCortexHeader(int childColor) {
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
        cc.setSampleName(GRAPH.getSampleName(childColor));

        ch.addColor(cc);

        return ch;
    }
}
