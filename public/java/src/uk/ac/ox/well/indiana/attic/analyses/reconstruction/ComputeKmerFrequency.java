package uk.ac.ox.well.indiana.attic.analyses.reconstruction;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.performance.PerformanceUtils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Map;

public class ComputeKmerFrequency extends Module {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public ArrayList<CortexGraph> CORTEX_GRAPHS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CortexKmer, Integer> kmerFreq = new java.util.HashMap<CortexKmer, Integer>();

        int highestFrequency = 0;

        for (CortexGraph cg : CORTEX_GRAPHS) {
            for (CortexRecord cr : cg) {
                CortexKmer ck = cr.getKmer();

                if (!kmerFreq.containsKey(ck)) {
                    kmerFreq.put(ck, 1);
                } else {
                    kmerFreq.put(ck, kmerFreq.get(ck) + 1);
                }

                if (kmerFreq.get(ck) > highestFrequency) {
                    highestFrequency = kmerFreq.get(ck);
                }
            }

            log.info("Loaded '{}'.  Map now contains {} kmers.  ({})", cg.getCortexFile().getAbsolutePath(), kmerFreq.size(), PerformanceUtils.getCompactMemoryUsageStats());
        }

        log.info("Computing frequency spectrum.");

        int[] freqCounts = new int[highestFrequency + 1];
        for (CortexKmer ck : kmerFreq.keySet()) {
            int freq = kmerFreq.get(ck);
            freqCounts[freq]++;
        }

        log.info("Writing kmer frequency to disk.");

        out.println("freq\tcount");
        for (int freq = 0; freq < freqCounts.length; freq++) {
            out.println(freq + "\t" + freqCounts[freq]);
        }
    }
}
