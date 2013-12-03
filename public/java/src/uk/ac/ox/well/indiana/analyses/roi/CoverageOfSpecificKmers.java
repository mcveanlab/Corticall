package uk.ac.ox.well.indiana.analyses.roi;

import com.google.common.base.Joiner;
import uk.ac.ox.well.indiana.tools.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

public class CoverageOfSpecificKmers extends Module {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="kmers", shortName="kmers", doc="Kmers to look for")
    public HashSet<String> KMERS;

    @Argument(fullName="normalizationGraph", shortName="n", doc="Normalize?")
    public CortexGraph NORMALIZATION_GRAPH;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, Map<String, Integer>> results = new TreeMap<String, Map<String, Integer>>();
        Map<String, Double> norms = new HashMap<String, Double>();

        // Find kmers that are unique in the genome
        Set<CortexRecord> unitCoverageKmers = new HashSet<CortexRecord>();
        for (CortexRecord cr : NORMALIZATION_GRAPH) {
            int cov = cr.getCoverage(0);

            if (cov == 1) {
                unitCoverageKmers.add(cr);
            }
        }

        int numRecords = 0;
        int index = 0;

        for (CortexRecord cr : CORTEX_GRAPH) {
            String fw = cr.getKmerAsString();
            String rc = SequenceUtils.reverseComplement(fw);
            String kmer = null;

            if (KMERS.contains(fw)) { kmer = fw; }
            else if (KMERS.contains(rc)) { kmer = rc; }

            if (kmer != null) {
                if (!results.containsKey(kmer)) {
                    results.put(kmer, new TreeMap<String, Integer>());
                }

                for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
                    String sample = CORTEX_GRAPH.getColor(color).getSampleName();

                    results.get(kmer).put(sample, cr.getCoverage(color));
                }
            }

            if (unitCoverageKmers.contains(cr.getKmer())) {
                for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
                    String sample = CORTEX_GRAPH.getColor(color).getSampleName();
                    int cov = cr.getCoverage(color);

                    if (!norms.containsKey(sample)) {
                        norms.put(sample, 0.0);
                    }

                    norms.put(sample, norms.get(sample) + cov);
                }

                numRecords++;
            }

            // Show our progress
            if (index % (CORTEX_GRAPH.getNumRecords() / 10) == 0) {
                log.info("processed {} records", index);
            }
            index++;
        }

        for (String sample : norms.keySet()) {
            double norm = norms.get(sample) / (double) numRecords;
            norms.put(sample, norm);
        }

        Set<String> samples = results.get(results.keySet().iterator().next()).keySet();

        out.println("\t" + Joiner.on("\t").join(samples));
        for (String kmer : results.keySet()) {
            List<String> fields = new ArrayList<String>();

            for (String sample : results.get(kmer).keySet()) {
                double coverage = results.get(kmer).get(sample);
                double normalizedCoverage = coverage / norms.get(sample);

                fields.add(String.valueOf(coverage) + ":" + String.valueOf(normalizedCoverage));
            }

            out.println(kmer + "\t" + Joiner.on("\t").join(fields));
        }
    }
}
