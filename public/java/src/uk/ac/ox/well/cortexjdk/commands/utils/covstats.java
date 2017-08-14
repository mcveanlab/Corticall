package uk.ac.ox.well.cortexjdk.commands.utils;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;

import java.io.PrintStream;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

public class covstats extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="child", shortName="c", doc="Child sample name")
    public String CHILD;

    @Argument(fullName="parent", shortName="p", doc="Parents")
    public HashSet<String> PARENTS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        Set<Integer> parentColors = getColorsForSampleNames(PARENTS);

        Map<Integer, Integer> hist = new TreeMap<>();
        long sharedRecords = 0L;

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing graph...")
                .message("records processed")
                .updateRecord(GRAPH.getNumRecords() / 10)
                .maxRecord(GRAPH.getNumRecords())
                .make(log);

        for (CortexRecord cr : GRAPH) {
            pm.update("records processed (" + sharedRecords + " shared records so far)");

            boolean isInChild = cr.getCoverage(childColor) > 0;
            int numberOfParents = 0;
            int numberOfChildren = 0;

            for (int c = 0; c < cr.getNumColors(); c++) {
                if (cr.getCoverage(c) > 0) {
                    if (childColor == c) { isInChild = true; }
                    else if (parentColors.contains(c)) { numberOfParents++; }
                    else { numberOfChildren++; }
                }
            }

            int childCov = cr.getCoverage(childColor);

            if (isInChild && numberOfParents > 0 && numberOfChildren > 0) {
                if (!hist.containsKey(childCov)) {
                    hist.put(childCov, 0);
                }

                hist.put(childCov, hist.get(childCov) + numberOfParents + numberOfChildren);

                sharedRecords++;
            }
        }

        for (int cov : hist.keySet()) {
            out.println(cov + "\t" + hist.get(cov));
        }
    }

    private Set<Integer> getColorsForSampleNames(Set<String> samples) {
        Set<Integer> colors = new HashSet<>();

        for (String sample : samples) {
            colors.add(GRAPH.getColorForSampleName(sample));
        }

        return colors;
    }
}
