package uk.ac.ox.well.cortexjdk.commands.prefilter;

import com.google.api.client.util.Joiner;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Description;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraphWriter;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

@Description(text="Remove kmers shared among children (as these are unlikely to tag de novo mutations)")
public class RemoveShared extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="ignore", shortName="i", doc="Ignore specified samples", required=false)
    public ArrayList<String> IGNORE;

    @Argument(fullName="roi", shortName="r", doc="ROIs")
    public CortexGraph ROI;

    @Output
    public File out;

    @Output(fullName="excluded_out", shortName="xo", doc="Excluded kmers output file")
    public File shared_out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        Set<Integer> parentColors = new HashSet<>(GRAPH.getColorsForSampleNames(PARENTS));
        Set<Integer> ignoreColors = new HashSet<>(GRAPH.getColorsForSampleNames(IGNORE));

        log.info("Colors:");
        log.info(" -   child: {}", GRAPH.getColorForSampleName(CHILD));
        log.info(" - parents: {}", GRAPH.getColorsForSampleNames(PARENTS));
        log.info(" -  ignore: {}", GRAPH.getColorsForSampleNames(IGNORE));

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Finding shared kmers")
                .message("records processed")
                .maxRecord(ROI.getNumRecords())
                .make(log);

        Set<CanonicalKmer> sharedKmers = new HashSet<>();

        for (CortexRecord rr : ROI) {
            if (!sharedKmers.contains(rr.getCanonicalKmer())) {
                CortexRecord cr = GRAPH.findRecord(rr.getCanonicalKmer());

                boolean isShared = false;

                for (int c = 0; c < GRAPH.getNumColors(); c++) {
                    if (c != childColor && !parentColors.contains(c) && !ignoreColors.contains(c) && cr.getCoverage(c) > 0) {
                        sharedKmers.add(rr.getCanonicalKmer());
                        isShared = true;

                        break;
                    }
                }

                if (isShared && log.isDebugEnabled()) {
                    List<String> records = new ArrayList<>();

                    for (int c = 0; c < GRAPH.getNumColors(); c++) {
                        records.add(String.format("%d:%d", c, cr.getCoverage(c)));
                    }

                    log.debug("{}", Joiner.on(' ').join(records));
                }
            }

            pm.update();
        }

        log.info("Found {} shared kmers", sharedKmers.size());

        log.info("Writing...");

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(ROI.getHeader());

        CortexGraphWriter cgo = new CortexGraphWriter(shared_out);
        cgo.setHeader(ROI.getHeader());

        int numKept = 0, numExcluded = 0;
        for (CortexRecord rr : ROI) {
            if (!sharedKmers.contains(rr.getCanonicalKmer())) {
                cgw.addRecord(rr);
                numKept++;
            } else {
                cgo.addRecord(rr);
                numExcluded++;
            }
        }

        cgw.close();
        cgo.close();

        log.info("  {}/{} ({}%) kept, {}/{} ({}%) excluded",
                numKept,     ROI.getNumRecords(), 100.0f * (float) numKept / (float) ROI.getNumRecords(),
                numExcluded, ROI.getNumRecords(), 100.0f * (float) numExcluded / (float) ROI.getNumRecords()
        );
    }
}
