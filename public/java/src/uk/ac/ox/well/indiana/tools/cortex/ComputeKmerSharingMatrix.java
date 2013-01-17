package uk.ac.ox.well.indiana.tools.cortex;

import com.google.common.base.Joiner;
import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.sketches.Sketch;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.ArrayList;

public class ComputeKmerSharingMatrix extends Tool {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="colors", shortName="c", doc="Colors to process")
    public ArrayList<Integer> COLORS;

    @Output
    public PrintStream out;

    public int execute() {
        ArrayList<String> header = new ArrayList<String>();
        header.add("");

        for (int color : COLORS) {
            header.add(CORTEX_GRAPH.getColor(color).getSampleName());
        }

        out.println(Joiner.on("\t").join(header));

        for (CortexRecord cr : CORTEX_GRAPH) {
            int[] coverages = cr.getCoverages();

            boolean hasCoverage = false;
            int numColorsWithKmer = 0;
            boolean hasZeroOrUnitCoverageInColors = true;

            for (int color : COLORS) {
                hasCoverage |= (coverages[color] > 0);
                numColorsWithKmer += coverages[color] == 0 ? 0 : 1;
                hasZeroOrUnitCoverageInColors &= (coverages[color] <= 1);
            }

            if (hasCoverage && numColorsWithKmer > 1 && hasZeroOrUnitCoverageInColors) {
                ArrayList<String> fields = new ArrayList<String>();
                fields.add(cr.getKmerString());

                for (int color : COLORS) {
                    fields.add(Integer.toString(coverages[color]));
                }

                out.println(Joiner.on("\t").join(fields));
            }
        }

        return 0;
    }
}
