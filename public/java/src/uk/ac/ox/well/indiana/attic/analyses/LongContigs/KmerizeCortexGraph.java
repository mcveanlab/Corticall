package uk.ac.ox.well.indiana.attic.analyses.LongContigs;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

public class KmerizeCortexGraph extends Module {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public ArrayList<CortexGraph> CORTEX_GRAPHS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Set<String> kmers = new HashSet<String>();

        log.info("Processing graphs...");

        for (CortexGraph cg : CORTEX_GRAPHS) {
            int counter = 0;
            for (CortexRecord cr : cg) {
                if (counter % (cg.getNumRecords()/10) == 0) {
                    log.info("  {}: processed {}/{} (~{}%) records", cg.getCortexFile().getName(), counter, cg.getNumRecords(), String.format("%.2f", 100.0*counter/cg.getNumRecords()));
                }
                counter++;

                String kmer = cr.getKmerAsString();
                kmers.add(kmer);
            }
        }

        log.info("Writing kmers to disk...");

        int counter = 1;
        for (String kmer : kmers) {
            out.println(">kmer." + counter);
            out.println(kmer);

            counter++;
        }
    }
}
