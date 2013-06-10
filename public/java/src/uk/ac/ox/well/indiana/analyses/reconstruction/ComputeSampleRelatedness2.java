package uk.ac.ox.well.indiana.analyses.reconstruction;

import com.google.common.base.Joiner;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataFrame;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.*;
import java.util.*;

public class ComputeSampleRelatedness2 extends Tool {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        DataFrame<String, String, Float> relatedness = new DataFrame<String, String, Float>(0.0f);

        int numRecords = 0;
        for (CortexRecord cr : CORTEX_GRAPH) {
            if (numRecords % (CORTEX_GRAPH.getNumRecords() / 10) == 0) {
                log.info("Processed {}/{} records", numRecords, CORTEX_GRAPH.getNumRecords());
            }
            numRecords++;

            for (int color1 = 0; color1 < CORTEX_GRAPH.getNumColors(); color1++) {
                for (int color2 = 0; color2 < CORTEX_GRAPH.getNumColors(); color2++) {
                    if (cr.getCoverage(color1) > 0 && cr.getCoverage(color2) > 0) {
                        String sample1 = CORTEX_GRAPH.getColor(color1).getSampleName();
                        String sample2 = CORTEX_GRAPH.getColor(color2).getSampleName();

                        float oldValue = relatedness.get(sample1, sample2);
                        relatedness.set(sample1, sample2, oldValue + 1.0f);
                    }
                }
            }
        }

        out.println(relatedness);
    }
}
