package uk.ac.ox.well.indiana.analyses.roi;

import uk.ac.ox.well.indiana.tools.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeMap;

public class CoverageOfSpecificKmers extends Module {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="kmers", shortName="kmers", doc="Kmers to look for")
    public HashSet<String> KMERS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, Map<String, Integer>> results = new TreeMap<String, Map<String, Integer>>();

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
        }

        out.println("kmer\tsample\tcoverage");
        for (String kmer : results.keySet()) {
            for (String sample : results.get(kmer).keySet()) {
                int coverage = results.get(kmer).get(sample);

                out.println(kmer + "\t" + sample + "\t" + coverage);
            }
        }
    }
}
