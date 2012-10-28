package uk.ac.ox.well.indiana.tools.cortex;

import net.sf.picard.reference.FastaSequenceFile;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;

import java.util.ArrayList;

public class ViewCortexCoverage extends ViewCortexBase {
    @Argument(fullName="targetsOfInterest", shortName="toi", doc="One or more fasta file of targets to find", required=false)
    public ArrayList<FastaSequenceFile> TARGETS_OF_INTEREST;

    @Override
    public int execute() {
//        log.info("Loading {}-mers from targets of interest", CORTEX_GRAPH.getKmerSize());
//        Map<Integer, String> kmerMap = null;
//        if (TARGETS_OF_INTEREST != null) {
//            kmerMap = SequenceUtils.loadSequenceCodesAsAlphanumericallyLowestKmers(TARGETS_OF_INTEREST, CORTEX_GRAPH.getKmerSize());
//        }

        log.info("Processing Cortex records");
        for (CortexRecord cr : CORTEX_GRAPH) {
            if (satisfiesConstraints(cr)) {
                out.println(cr);
            }
        }

        return 0;
    }
}
