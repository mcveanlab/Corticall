package uk.ac.ox.well.indiana.analyses.reconstruction;

import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;

import java.io.PrintStream;
import java.util.ArrayList;

public class ExtractCortexSampleNames extends Tool {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public ArrayList<CortexGraph> CORTEX_GRAPHS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        for (CortexGraph cg : CORTEX_GRAPHS) {
            for (int color = 0; color < cg.getNumColors(); color++) {
                String sample = cg.getColor(color).getSampleName();

                out.println(sample);
            }
        }
    }
}
