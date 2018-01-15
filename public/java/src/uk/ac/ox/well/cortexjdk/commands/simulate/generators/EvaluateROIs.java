package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

public class EvaluateROIs extends Module {
    @Argument(fullName="knownRois", shortName="k", doc="Known ROIs")
    public File KNOWN_ROIS;

    @Argument(fullName="rois", shortName="r", doc="ROIs")
    public CortexGraph ROIS;

    @Override
    public void execute() {
        TableReader tr = new TableReader(KNOWN_ROIS, "id", "numKmers", "index", "kmer");

        Map<CanonicalKmer, Map<String, String>> knownRois = new HashMap<>();
        Map<CanonicalKmer, Boolean> found = new LinkedHashMap<>();
        for (Map<String, String> te : tr) {
            CanonicalKmer ck = new CanonicalKmer(te.get("kmer"));
            knownRois.put(ck, te);
            found.put(ck, false);
        }

        for (CortexRecord cr : ROIS) {
            if (found.containsKey(cr.getCanonicalKmer())) {
                found.put(cr.getCanonicalKmer(), true);
            } else {
                log.info("missing: {}", cr);
            }
        }

        for (CanonicalKmer ck : found.keySet()) {
            Map<String, String> te = knownRois.get(ck);

            log.info("{} {} {} {} {}", te.get("id"), te.get("numKmers"), te.get("index"), te.get("kmer"), found.get(ck));
        }
    }
}
