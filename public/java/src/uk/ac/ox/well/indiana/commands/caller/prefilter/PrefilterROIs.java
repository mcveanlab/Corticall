package uk.ac.ox.well.indiana.commands.caller.prefilter;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

@Description(text="Remove obvious sequencing errors")
public class PrefilterROIs extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="roi", shortName="r", doc="ROIs")
    public CortexGraph ROI;

    @Output
    public File out;

    @Override
    public void execute() {
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(PARENTS);
        int childColor = GRAPH.getColorForSampleName(CHILD);

        Set<Integer> usedColors = new HashSet<>();
        usedColors.addAll(parentColors);
        usedColors.add(childColor);

        int numOtherChildren = GRAPH.getNumColors() - usedColors.size();

        for (CortexRecord rr : ROI) {
            CortexRecord cr = GRAPH.findRecord(rr.getCortexKmer());

            if (cr == null) {
                throw new IndianaException("Found an ROI record that's not in the mDBG (" + rr + ").");
            }

            int presentInChildren = 0;

            for (int c = 0; c < GRAPH.getNumColors(); c++) {
                if (!usedColors.contains(c) && cr.getCoverage(c) > 0) {
                    presentInChildren++;
                }
            }


        }
    }
}
