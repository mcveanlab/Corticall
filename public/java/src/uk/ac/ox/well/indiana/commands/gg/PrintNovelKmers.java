package uk.ac.ox.well.indiana.commands.gg;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexColor;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

public class PrintNovelKmers extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="samplesToConsider", shortName="s", doc="Samples to consider (regex)")
    public ArrayList<String> SAMPLE_REGEX;

    @Output
    public PrintStream out;

    private boolean isNovel(CortexGraph cg, CortexRecord cr, Set<String> samplesToConsider) {
        boolean consideredSamplesHaveCoverage = false;
        boolean otherSamplesHaveCoverage = false;

        for (int c = 0; c < cg.getNumColors(); c++) {
            CortexColor cc = cg.getColor(c);

            if (cr.getCoverage(c) > 0) {
                if (samplesToConsider.contains(cc.getSampleName())) {
                    consideredSamplesHaveCoverage = true;
                } else {
                    otherSamplesHaveCoverage = true;
                }
            }
        }

        return consideredSamplesHaveCoverage && !otherSamplesHaveCoverage;
    }

    @Override
    public void execute() {
        Set<String> samplesToConsider = new HashSet<String>();
        for (String regex : SAMPLE_REGEX) {
            for (CortexColor cc : GRAPH.getColors()) {
                if (cc.getSampleName().matches(regex)) {
                    samplesToConsider.add(cc.getSampleName());
                }
            }
        }

        log.info("Considering novel kmers in the following samples:");
        for (String sample : samplesToConsider) {
            log.info("  {}", sample);
        }

        log.info("Looking for novel kmers...");
        int numKmers = 0, numNovelKmers = 0;
        for (CortexRecord cr : GRAPH) {
            if (numKmers % (GRAPH.getNumRecords() / 10) == 0) {
                log.info("  {}/{} kmers processed", numKmers, GRAPH.getNumRecords());
            }
            numKmers++;

            if (isNovel(GRAPH, cr, samplesToConsider)) {
                numNovelKmers++;
            }
        }

        log.info("  found {} novel kmers", numNovelKmers);
    }
}
