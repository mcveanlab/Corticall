package uk.ac.ox.well.indiana.attic.analyses.geo;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.assembly.cortex.CortexGraphWalker;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexMap;

import java.io.PrintStream;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class GetContigs extends Module {
    @Argument(fullName="diagnosticKmers", shortName="dk", doc="Diagnostic kmers")
    public HashSet<CortexKmer> DIAGNOSTIC_KMERS;

    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexMap CORTEX_MAP;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        CortexGraphWalker cgw = new CortexGraphWalker(CORTEX_MAP);

        for (int color = 0; color < CORTEX_MAP.getGraph().getNumColors(); color++) {
            String sampleName = CORTEX_MAP.getGraph().getColor(color).getSampleName();

            Map<CortexKmer, Set<CortexKmer>> contigs = cgw.buildContigs(color, DIAGNOSTIC_KMERS);

            int index = 1;
            for (CortexKmer contig : contigs.keySet()) {
                out.println(">" + sampleName + "." + index);
                out.println(contig);

                index++;
            }
        }
    }
}
