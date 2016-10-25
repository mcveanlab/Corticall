package uk.ac.ox.well.indiana.commands.playground.kmerdb;

import com.google.common.base.Joiner;
import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.collection.CortexCollection2;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.*;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.File;
import java.util.*;

public class EmitDbNovelKeys extends Module {
    @Argument(fullName="graphs", shortName="g", doc="Graphs")
    public CortexCollection2 GRAPHS;

    @Argument(fullName="childColor", shortName="cc", doc="Child color")
    public Integer CHILD_COLOR = 0;

    @Argument(fullName="parentColor", shortName="pc", doc="Parent color")
    public HashSet<Integer> PARENT_COLORS = new HashSet<>();

    @Argument(fullName="coverageLowerLimit", shortName="l", doc="Coverage lower limit")
    public Integer COVERAGE_LOWER_LIMIT = 0;

    @Output
    public File out;

    @Override
    public void execute() {
        if (PARENT_COLORS.isEmpty()) {
            for (int c = 0; c < GRAPHS.getNumColors(); c++) {
                if (c != CHILD_COLOR) {
                    PARENT_COLORS.add(c);
                }
            }
        }

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(makeCortexHeader());

        log.info("Colors: (child) {} (parents) {}", CHILD_COLOR, Joiner.on(",").join(PARENT_COLORS));

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing graph...")
                .message("records processed")
                .updateRecord(GRAPHS.getGraph(0).getNumRecords() / 10)
                .make(log);

        long numNovelRecords = 0;
        for (CortexRecord cr : GRAPHS) {
            if (CortexUtils.isNovelKmer(cr, CHILD_COLOR, PARENT_COLORS) && cr.getCoverage(CHILD_COLOR) > COVERAGE_LOWER_LIMIT) {
                CortexRecord novelCr = new CortexRecord(
                        cr.getBinaryKmer(),
                        new int[] { cr.getCoverages()[CHILD_COLOR] },
                        new byte[] { cr.getEdges()[CHILD_COLOR] },
                        cr.getKmerSize(),
                        cr.getKmerBits()
                );

                cgw.addRecord(novelCr);

                numNovelRecords++;
            }

            pm.update("records processed (" + numNovelRecords + " novel so far)");
        }

        log.info("  {} total novel records", numNovelRecords);

        cgw.close();
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
        cc.setSampleName("novel");

        ch.addColor(cc);
        return ch;
    }
}
