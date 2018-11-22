package uk.ac.ox.well.cortexjdk.commands.prefilter;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraphWriter;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

public class FindLowComplexity extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="crThreshold", shortName="t", doc="Complexity threshold")
    public Float COMPLEXITY_THRESHOLD = 0.70f;

    @Output
    public File out;

    @Override
    public void execute() {
        String child = ROI.getSampleName(0);

        //int childColor = GRAPH.getColorForSampleName(child);
        //Set<Integer> parentColors = new HashSet<>(GRAPH.getColorsForSampleNames(PARENTS));

        log.info("Colors:");
        log.info(" -   child: {}", GRAPH.getColorForSampleName(child));
        log.info(" - parents: {}", GRAPH.getColorsForSampleNames(PARENTS));

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Finding low complexity kmers")
                .message("records processed")
                .maxRecord(ROI.getNumRecords())
                .make(log);

        Set<CanonicalKmer> lowComplexity = new HashSet<>();

        for (CortexRecord rr : ROI) {
            //log.info("{} {} {}", rr.getCanonicalKmer(), SequenceUtils.computeCompressionRatio(rr.getCanonicalKmer()), isLowComplexity(rr, COMPLEXITY_THRESHOLD));

            if (isLowComplexity(rr, COMPLEXITY_THRESHOLD)) {
                lowComplexity.add(rr.getCanonicalKmer());

            }

            pm.update();
        }

        log.info("Found {} low complexity kmers", lowComplexity.size());

        log.info("Writing...");

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(ROI.getHeader());

        //CortexGraphWriter cgo = new CortexGraphWriter(dust_out);
        //cgo.setHeader(ROI.getHeader());

        int numKept = 0, numExcluded = 0;
        for (CortexRecord rr : ROI) {
            if (!lowComplexity.contains(rr.getCanonicalKmer())) {
                //cgw.addRecord(rr);
                numKept++;
            } else {
                //cgo.addRecord(rr);
                cgw.addRecord(rr);
                numExcluded++;
            }
        }

        cgw.close();
        //cgo.close();

        log.info("  {}/{} ({}%) kept, {}/{} ({}%) excluded",
                numKept,     ROI.getNumRecords(), 100.0f * (float) numKept / (float) ROI.getNumRecords(),
                numExcluded, ROI.getNumRecords(), 100.0f * (float) numExcluded / (float) ROI.getNumRecords()
        );
    }

    private boolean isLowComplexity(CortexRecord cr, float cratioThreshold) {
        float cratio = SequenceUtils.computeCompressionRatio(cr.getCanonicalKmer());

        return cratio < cratioThreshold;
    }
}
