package uk.ac.ox.well.indiana.analyses.kmerSharing;

import com.google.common.base.Joiner;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.assembly.CortexGraphWalker;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.utils.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class GetSupernodes extends Tool {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="kmerReferencePanel", shortName="krp", doc="Kmer reference panel")
    public File KMER_REFERENCE_PANEL;

    @Output
    public PrintStream out;

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

    private Map<String, String> loadKmerReferencePanel() {
        Map<String, String> panel = new HashMap<String, String>();

        TableReader tr = new TableReader(KMER_REFERENCE_PANEL);
        for (HashMap<String, String> entry : tr) {
            panel.put(entry.get("kmer"), entry.get("genes"));
        }

        if (panel.containsKey("GCACGAAGTTTTGCAGATATAGGCGATATTA")) {
            log.info("Panel contains kmer: {}", panel.get("GCACGAAGTTTTGCAGATATAGGCGATATTA"));
        } else {
            log.info("Panel is missing kmer???");
        }

        return panel;
    }

    @Override
    public int execute() {
        Map<String, String> panel = loadKmerReferencePanel();
        Map<String, CortexRecord> records = loadCortexRecords();

        CortexGraphWalker cgw = new CortexGraphWalker(records);

        String[] header = { "color", "genes", "superNode" };
        out.println(Joiner.on("\t").join(header));

        boolean debug = false;

        for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
            Collection<String> supernodes = cgw.getReferenceGuidedSupernodes(color, panel.keySet());

            for (String supernode : supernodes) {
                if (supernode.equalsIgnoreCase("GCACGAAGTTTTGCAGATATAGGCGATATTA")) {
                    debug = true;
                }

                Set<String> genes = new HashSet<String>();

                for (int i = 0; i <= supernode.length() - CORTEX_GRAPH.getKmerSize(); i++) {
                    String fw = SequenceUtils.getAlphanumericallyLowestOrientation(supernode.substring(i, i + CORTEX_GRAPH.getKmerSize()));

                    if (panel.containsKey(fw)) {
                        genes.add(panel.get(fw));

                        if (debug) {
                            log.debug("fw: {} {}", fw, panel.get(fw));
                        }
                    } else {
                        if (debug) {
                            log.debug("fw: {} missing", fw);
                        }
                    }
                }

                String[] entry = { String.valueOf(color), Joiner.on(",").join(genes), supernode };
                out.println(Joiner.on("\t").join(entry));

                debug = false;
            }
        }

        return 0;
    }
}
