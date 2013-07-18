package uk.ac.ox.well.indiana.analyses.kmerSharing;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import uk.ac.ox.well.indiana.tools.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

public class FindZeroThresholds extends Module {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public File CORTEX_GRAPH;

    @Argument(fullName="reference", shortName="R", doc="Reference file")
    public IndexedFastaSequenceFile REFERENCE;

    @Argument(fullName="gff", shortName="gff", doc="GFF file")
    public GFF3 GFF;

    private Map<String, String> loadSequences() {
        Map<String, String> sequences = new HashMap<String, String>();

        for (GFF3Record r : GFF) {
            if ("gene".equalsIgnoreCase(r.getType())) {
                String seq = new String(REFERENCE.getSubsequenceAt(r.getSeqid(), r.getStart(), r.getEnd()).getBases());

                sequences.put(r.getAttribute("ID"), seq);
            }
        }

        return sequences;
    }

    @Override
    public void execute() {
        Map<String, String> sequences = loadSequences();

        for (String geneName1 : sequences.keySet()) {
            for (String geneName2 : sequences.keySet()) {

            }
        }
    }
}
