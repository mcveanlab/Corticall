package uk.ac.ox.well.cortexjdk.utils.samples;

import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;
import uk.ac.ox.well.cortexjdk.utils.io.utils.LineReader;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created by kiran on 16/08/2017.
 */
public class SampleOrganizer {
    private File pedFile;
    private Map<String, SampleInfo> siMap = new HashMap<>();

    public SampleOrganizer(String ped) {
        loadPed(new File(ped));
    }

    public SampleOrganizer(File pedFile) {
        loadPed(pedFile);
    }

    private void loadPed(File pedFile) {
        this.pedFile = pedFile;

        TableReader tr = new TableReader(pedFile, "family", "sample", "mother", "father", "NA1", "NA2");

        for (Map<String, String> te : tr) {
            String sample = te.get("sample");

            String mother = te.get("mother");
            if (mother.equals("0")) { mother = null; }

            String father = te.get("father");
            if (father.equals("0")) { father = null; }

            if ((mother != null || father != null) && !siMap.containsKey(sample)) {
                siMap.put(sample, new SampleInfo(sample, mother, father));
            }
        }

        for (Map<String, String> te : tr) {
            String sample = te.get("sample");

            String mother = te.get("mother");
            if (mother.equals("0")) { mother = null; }

            String father = te.get("father");
            if (father.equals("0")) { father = null; }

            if (mother == null && father == null) {
                for (SampleInfo si : siMap.values()) {

                }
            }
        }

    }
}
