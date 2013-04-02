package uk.ac.ox.well.indiana.analyses.kmerSharing;

import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.assembly.CortexGraphWalker;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.utils.TableReader;

import java.io.File;
import java.util.*;

public class GetSupernodes extends Tool {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="kmerReferencePanel", shortName="krp", doc="Kmer reference panel")
    public File KMER_REFERENCE_PANEL;

    private Map<String, CortexRecord> loadCortexRecords() {
        Map<String, CortexRecord> records = new HashMap<String, CortexRecord>();

        int recordNum = 0;
        for (CortexRecord cr : CORTEX_GRAPH) {
            if (recordNum % (CORTEX_GRAPH.getNumRecords() / 10) == 0) {
                log.info("Loaded {}/{} records", recordNum, CORTEX_GRAPH.getNumRecords());
            }
            recordNum++;

            String kmer = cr.getKmerString();

            records.put(kmer, cr);
        }

        return records;
    }

    private Set<String> loadKmerReferencePanel() {
        Set<String> panel = new HashSet<String>();

        TableReader tr = new TableReader(KMER_REFERENCE_PANEL);
        for (HashMap<String, String> entry : tr) {
            panel.add(entry.get("kmer"));
        }

        return panel;
    }

    @Override
    public int execute() {
        Map<String, CortexRecord> records = loadCortexRecords();
        Set<String> panel = loadKmerReferencePanel();

        CortexGraphWalker cgw = new CortexGraphWalker(records);
        Collection<String> supernodes = cgw.getReferenceGuidedSupernodes(0, panel);

        log.info("Found {} supernodes", supernodes.size());
        for (String supernode : supernodes) {
            log.info("{} {}", supernode.length(), supernode);
        }

        return 0;
    }
}
