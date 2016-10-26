package uk.ac.ox.well.indiana.commands.playground.kmerdb;

import com.google.common.base.Joiner;
import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.collection.CortexCollection;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.*;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.File;
import java.util.*;

public class EmitDbNovelKeys extends Module {
    @Argument(fullName="graphs", shortName="g", doc="Graphs")
    public CortexCollection GRAPHS;

    @Argument(fullName="childColor", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="parentColor", shortName="p", doc="Parent")
    public HashSet<String> PARENTS = new HashSet<>();

    @Argument(fullName="coverageLowerLimit", shortName="l", doc="Coverage lower limit")
    public Integer COVERAGE_LOWER_LIMIT = 0;

    @Output
    public File out;

    @Override
    public void execute() {
        log.info("Colors:");

        int child_color = 0;
        Set<Integer> parent_colors = new TreeSet<>();

        for (int c = 0; c < GRAPHS.getNumColors(); c++) {
            if (CHILD.contains(GRAPHS.getSampleName(c))) {
                child_color = c;

                log.info("  child: {} {}", CHILD, child_color);

                break;
            }
        }

        for (int c = 0; c < GRAPHS.getNumColors(); c++) {
            if (PARENTS.contains(GRAPHS.getSampleName(c))) {
                parent_colors.add(c);

                log.info("  parent: {} {}", GRAPHS.getSampleName(c), c);
            }
        }

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(makeCortexHeader());

        //log.info("Colors: (child) {} (parents) {}", child_color, Joiner.on(",").join(parent_colors));

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing graph...")
                .message("records processed")
                .updateRecord(GRAPHS.getGraph(0).getNumRecords() / 10)
                .make(log);

        long numNovelRecords = 0;
        for (CortexRecord cr : GRAPHS) {
            if (CortexUtils.isNovelKmer(cr, child_color, parent_colors) && cr.getCoverage(child_color) > COVERAGE_LOWER_LIMIT) {
                CortexRecord novelCr = new CortexRecord(
                        cr.getBinaryKmer(),
                        new int[] { cr.getCoverages()[child_color] },
                        new byte[] { cr.getEdges()[child_color] },
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
