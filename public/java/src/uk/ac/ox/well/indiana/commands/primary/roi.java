package uk.ac.ox.well.indiana.commands.primary;

import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.*;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;

import java.io.File;

@Description(text="Identify regions of interest (putative de novo mutations) in graphs")
public class roi extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPHS;

    @Argument(fullName="childColor", shortName="c", doc="Child")
    public String CHILD;

    @Output
    public File out;

    @Override
    public void execute() {
        int childColor = GRAPHS.getColorForSampleName(CHILD);

        log.info("Color: {} {}", CHILD, childColor);

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing graph...")
                .message("records processed")
                .updateRecord(GRAPHS.getNumRecords() / 10)
                .maxRecord(GRAPHS.getNumRecords())
                .make(log);

        long numNovelRecords = 0L;

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(makeCortexHeader());

        for (CortexRecord cr : GRAPHS) {
            if (isNovel(cr, childColor)) {
                CortexRecord novelCr = new CortexRecord(
                    cr.getBinaryKmer(),
                    new int[] { cr.getCoverages()[childColor] },
                    new byte[] { cr.getEdges()[childColor] },
                    cr.getKmerSize(),
                    cr.getKmerBits()
                );

                cgw.addRecord(novelCr);

                numNovelRecords++;
            }

            pm.update("records processed (" + numNovelRecords + " novel so far)");
        }

        cgw.close();
    }

    private boolean isNovel(CortexRecord cr, int childColor) {
        boolean childHasCoverage = cr.getCoverage(childColor) > 0;
        boolean othersHaveCoverage = false;

        for (int c = 0; c < cr.getNumColors(); c++) {
            if (c != childColor && cr.getCoverage(c) > 0) {
                othersHaveCoverage = true;
                break;
            }
        }

        return childHasCoverage && !othersHaveCoverage;
    }

    @NotNull
    private CortexHeader makeCortexHeader() {
        CortexHeader ch = new CortexHeader();
        ch.setVersion(6);
        ch.setNumColors(1);
        ch.setKmerSize(GRAPHS.getKmerSize());
        ch.setKmerBits(GRAPHS.getKmerBits());

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
