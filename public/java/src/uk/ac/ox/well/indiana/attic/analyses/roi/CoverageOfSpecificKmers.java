package uk.ac.ox.well.indiana.attic.analyses.roi;

import com.google.common.base.Joiner;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

public class CoverageOfSpecificKmers extends Module {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="normalizationGraph", shortName="n", doc="Normalize?")
    public CortexGraph NORMALIZATION_GRAPH;

    @Argument(fullName="kmers", shortName="kmers", doc="Kmers to look for")
    public HashSet<String> KMERS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, Map<String, Integer>> results = new TreeMap<String, Map<String, Integer>>();
        Map<String, Double> norms = new HashMap<String, Double>();

        // Find kmers that are unique in the genome
        log.info("Loading kmers with copy-number 1...");

        Set<CortexKmer> unitCoverageKmers = new HashSet<CortexKmer>();
        for (CortexRecord cr : NORMALIZATION_GRAPH) {
            int cov = cr.getCoverage(0);

            if (cov == 1) {
                unitCoverageKmers.add(cr.getKmer());
            }
        }

        log.info("Found {} unique kmers in the genome", unitCoverageKmers.size());

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
            if (index % (CORTEX_GRAPH.getNumRecords() / 5) == 0) {
                log.info("processed {}/{} records", index, CORTEX_GRAPH.getNumRecords());
            }
            index++;
        }

        for (String sample : norms.keySet()) {
            double norm = norms.get(sample) / (double) numRecords;

            log.info("normalization: sample={} cov={} records={} norm={}", sample, norms.get(sample), numRecords, norm);

            norms.put(sample, norm);
        }

        Set<String> samples = results.get(results.keySet().iterator().next()).keySet();
        Map<String, Double> totals = new HashMap<String, Double>();

        out.println("\t" + Joiner.on("\t").join(samples));
        for (String kmer : results.keySet()) {
            List<String> fields = new ArrayList<String>();

            for (String sample : samples) {
                double coverage = results.get(kmer).get(sample);
                double normalizedCoverage = coverage / norms.get(sample);

                //fields.add(String.valueOf(coverage) + ":" + String.valueOf(normalizedCoverage));
                fields.add(String.valueOf(normalizedCoverage));

                if (!totals.containsKey(sample)) { totals.put(sample, 0.0); }
                totals.put(sample, totals.get(sample) + normalizedCoverage);
            }

            out.println(kmer + "\t" + Joiner.on("\t").join(fields));
        }

        List<Double> fields = new ArrayList<Double>();
        for (String sample : samples) {
            double total = totals.get(sample);

            fields.add(total);
        }

        out.println("total" + "\t" + Joiner.on("\t").join(fields));
    }
}
