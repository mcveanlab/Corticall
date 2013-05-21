package uk.ac.ox.well.indiana.analyses.reconstruction;

import com.google.common.base.Joiner;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

/**
 * For each gene, show how many unique kmers were recovered per sample
 */
public class ComputeKmerRecoveryPerGene extends Tool {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="kmerReferencePanel", shortName="krp", doc="Kmer reference panel")
    public File KMER_REFERENCE_PANEL;

    @Output
    public PrintStream out;

    private class KmerObs {
        //         sample  counts
        public Map<String, Integer> seen = new TreeMap<String, Integer>();
        public Integer total = 0;
    }

    private Map<CortexKmer, String> loadKmerReferencePanel() {
        Map<CortexKmer, String> kmerReferencePanel = new HashMap<CortexKmer, String>();

        TableReader tr = new TableReader(KMER_REFERENCE_PANEL);
        for (Map<String, String> te : tr) {
            CortexKmer kmer = new CortexKmer(te.get("kmer"));
            String gene = te.get("gene");

            kmerReferencePanel.put(kmer, gene);
        }

        return kmerReferencePanel;
    }

    @Override
    public void execute() {
        //  kmers       gene
        Map<CortexKmer, String> krp = loadKmerReferencePanel();

        //  kmer        samples
        Map<CortexKmer, CortexRecord> kmerPresence = new HashMap<CortexKmer, CortexRecord>();

        int recordNum = 0;
        int storedRecords = 0;
        for (CortexRecord cr : CORTEX_GRAPH) {
            if (recordNum % (CORTEX_GRAPH.getNumRecords() / 10) == 0) {
                log.info("Processed {}/{} records ({} stored)", recordNum, CORTEX_GRAPH.getNumRecords(), storedRecords);
            }
            recordNum++;

            CortexKmer ck = cr.getKmer();

            if (krp.containsKey(ck)) {
                kmerPresence.put(ck, cr);
                storedRecords++;
            }
        }

        //  gene    KmerObs
        Map<String, KmerObs> kmerObs = new TreeMap<String, KmerObs>();

        for (CortexKmer ck : krp.keySet()) {
            String geneName = krp.get(ck);

            if (!kmerObs.containsKey(geneName)) {
                kmerObs.put(geneName, new KmerObs());
            }

            if (kmerPresence.containsKey(ck)) {
                kmerObs.get(geneName).total++;

                for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
                    String sampleName = CORTEX_GRAPH.getColor(color).getSampleName();

                    if (!kmerObs.get(geneName).seen.containsKey(sampleName)) {
                        kmerObs.get(geneName).seen.put(sampleName, 0);
                    }

                    CortexRecord cr = kmerPresence.get(ck);

                    if (cr.getCoverage(color) > 0) {
                        int oldCount = kmerObs.get(geneName).seen.get(sampleName);

                        kmerObs.get(geneName).seen.put(sampleName, oldCount + 1);
                    }
                }
            }
        }

        List<String> header = new ArrayList<String>();
        header.add("gene");
        header.add("total");

        for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
            header.add(CORTEX_GRAPH.getColor(color).getSampleName());
        }

        out.println(Joiner.on("\t").join(header));

        for (String geneName : kmerObs.keySet()) {
            List<String> fields = new ArrayList<String>();

            fields.add(geneName);
            fields.add(String.valueOf(kmerObs.get(geneName).total));

            for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
                String sampleName = CORTEX_GRAPH.getColor(color).getSampleName();

                int count = 0;
                if (kmerObs.get(geneName).seen.containsKey(sampleName)) {
                    count = kmerObs.get(geneName).seen.get(sampleName);
                }

                fields.add(String.valueOf(count));
            }

            out.println(Joiner.on("\t").join(fields));
        }
    }
}
