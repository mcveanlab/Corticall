package uk.ac.ox.well.cortexjdk.commands.simulate;

import com.google.common.base.Joiner;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class EvaluateROIs extends Module {
    @Argument(fullName="knownRois", shortName="k", doc="Known ROIs")
    public File KNOWN_ROIS;

    @Argument(fullName="rois", shortName="r", doc="ROIs")
    public CortexGraph ROIS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Set<CanonicalKmer> allRois = new HashSet<>();

        TableReader tr = new TableReader(KNOWN_ROIS, "id", "numKmers", "index", "kmer");

        Map<CanonicalKmer, Map<String, String>> knownRois = new HashMap<>();
        Map<CanonicalKmer, Boolean> found = new LinkedHashMap<>();
        for (Map<String, String> te : tr) {
            CanonicalKmer ck = new CanonicalKmer(te.get("kmer"));
            knownRois.put(ck, te);
            found.put(ck, false);

            allRois.add(ck);
        }

        Set<CanonicalKmer> calledRois = new HashSet<>();
        for (CortexRecord cr : ROIS) {
            allRois.add(cr.getCanonicalKmer());
            calledRois.add(cr.getCanonicalKmer());

            if (found.containsKey(cr.getCanonicalKmer())) {
                found.put(cr.getCanonicalKmer(), true);
            }
        }

        int uniqueToRois = 0, shared = 0, uniqueToKnowns = 0;

        out.println(Joiner.on("\t").join("class", "ck", "id", "numKmers", "index", "kmer", "found"));

        for (CanonicalKmer ck : allRois) {
            if (!knownRois.containsKey(ck)) {
                log.info("uniqueToRois   {}", ck);
                uniqueToRois++;

                out.println(Joiner.on("\t").join("uniqueToRois", ck, "NA", "NA", "NA", "NA", "NA"));
            } else if (found.containsKey(ck) && found.get(ck)) {
                Map<String, String> te = knownRois.get(ck);

                log.info("shared         {} {} {} {} {} {}", ck, te.get("id"), te.get("numKmers"), te.get("index"), te.get("kmer"), found.get(ck));
                shared++;

                out.println(Joiner.on("\t").join("shared", ck, te.get("id"), te.get("numKmers"), te.get("index"), te.get("kmer"), found.get(ck)));
            } else {
                Map<String, String> te = knownRois.get(ck);

                log.info("uniqueToKnowns {} {} {} {} {} {}", ck, te.get("id"), te.get("numKmers"), te.get("index"), te.get("kmer"), found.get(ck));
                uniqueToKnowns++;

                out.println(Joiner.on("\t").join("uniqueToKnowns", ck, te.get("id"), te.get("numKmers"), te.get("index"), te.get("kmer"), found.get(ck)));
            }
        }

        log.info("uniqueToRois={}, shared={}, uniqueToKnowns={}", uniqueToRois, shared, uniqueToKnowns);
    }
}
