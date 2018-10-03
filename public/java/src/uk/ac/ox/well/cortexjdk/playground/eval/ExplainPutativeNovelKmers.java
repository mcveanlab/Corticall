package uk.ac.ox.well.cortexjdk.playground.eval;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

public class ExplainPutativeNovelKmers extends Module {
    @Argument(fullName="rois", shortName="r", doc="ROIs")
    public CortexGraph GRAPH;

    @Argument(fullName="bam", shortName="b", doc="BAM")
    public SAMFileReader BAM;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CanonicalKmer, Boolean> rois = new HashMap<>();

        for (CortexRecord cr : GRAPH) {
            if (cr.getCoverage(0) == 0 && cr.getCoverage(1) == 0 && cr.getCoverage(2) > 0) {
                rois.put(cr.getCanonicalKmer(), null);
            }
        }

        for (SAMRecord sr : BAM) {
            for (int i = 0; i <= sr.getReadString().length() - GRAPH.getKmerSize(); i++) {
                CanonicalKmer ck = new CanonicalKmer(sr.getReadString().substring(i, i + GRAPH.getKmerSize()));

                if (rois.containsKey(ck) && rois.get(ck) == null) {
                    rois.put(ck, sr.getMappingQuality() > 0);
                }
            }
        }

        int mappedUniquely = 0;
        int mappedNonuniquely = 0;

        for (CanonicalKmer ck : rois.keySet()) {
            if (rois.get(ck) != null) {
                if (rois.get(ck)) {
                    mappedUniquely++;
                } else {
                    mappedNonuniquely++;
                }
            }
        }

        out.println(Joiner.on("\t").join(GRAPH.getSampleName(2), mappedUniquely, mappedNonuniquely));
    }
}
