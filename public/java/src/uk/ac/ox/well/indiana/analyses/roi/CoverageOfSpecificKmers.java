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

    @Argument(fullName="normalize", shortName="n", doc="Normalize?")
    public Boolean NORMALIZE = false;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, Map<String, Integer>> results = new TreeMap<String, Map<String, Integer>>();
        Map<String, Long> coverages = new HashMap<String, Long>();

        int numRecords = 0;
        int index = 0;

        for (CortexRecord cr : CORTEX_GRAPH) {
            String fw = cr.getKmerAsString();
            String rc = SequenceUtils.reverseComplement(fw);

            String kmer = null;

            if (KMERS.contains(fw)) {
                kmer = fw;
            } else if (KMERS.contains(rc)) {
                kmer = rc;
            }

            if (kmer != null) {
                if (!results.containsKey(kmer)) {
                    results.put(kmer, new TreeMap<String, Integer>());
                }

                for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
                    String sample = CORTEX_GRAPH.getColor(color).getSampleName();

                    results.get(kmer).put(sample, cr.getCoverage(color));
                }
            }

            if (index % 100 == 0) {
                for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
                    String sample = CORTEX_GRAPH.getColor(color).getSampleName();
                    int cov = cr.getCoverage(color);

                    if (!coverages.containsKey(sample)) {
                        coverages.put(sample, 0l);
                    }

                    coverages.put(sample, coverages.get(sample) + cov);
                }

                numRecords++;
            }

            index++;
        }

        Set<String> samples = results.get(results.keySet().iterator().next()).keySet();

        out.println("\t" + Joiner.on("\t").join(samples));
        for (String kmer : results.keySet()) {
            List<Float> fields = new ArrayList<Float>();

            for (String sample : results.get(kmer).keySet()) {
                float coverage = results.get(kmer).get(sample);
                float average = (float) coverages.get(sample) / (float) numRecords;
                float norm = (NORMALIZE) ? average : coverage;

                fields.add(norm);
            }

            out.println(kmer + "\t" + Joiner.on("\t").join(fields));
        }
    }
}
