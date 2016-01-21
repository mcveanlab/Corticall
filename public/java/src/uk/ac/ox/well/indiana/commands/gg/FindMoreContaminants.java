package uk.ac.ox.well.indiana.commands.gg;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraphWriter;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Set;

public class FindMoreContaminants extends Module {
    @Argument(fullName="clean", shortName="c", doc="Graph")
    public CortexGraph CLEAN;

    @Argument(fullName="contaminants", shortName="x", doc="Contaminants")
    public CortexGraph REJECTED_KMERS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Set<CortexKmer> contaminatingKmers = new HashSet<CortexKmer>();

        log.info("Exploring contaminants...");
        for (CortexRecord cr : REJECTED_KMERS) {
            log.info("{}", cr);

            if (!contaminatingKmers.contains(cr.getCortexKmer())) {
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs = CortexUtils.dfs(CLEAN, null, cr.getKmerAsString(), 0, null, ContaminantStopper.class);

                for (AnnotatedVertex rv : dfs.vertexSet()) {
                    CortexRecord rr = CLEAN.findRecord(new CortexKmer(rv.getKmer()));

                    if (rr.getCoverage(0) > 0 && rr.getCoverage(1) == 0 && rr.getCoverage(2) == 0) {
                        contaminatingKmers.add(rr.getCortexKmer());
                    }
                }

                contaminatingKmers.add(cr.getCortexKmer());
            }
        }

        int index = 0;
        for (CortexKmer ck : contaminatingKmers) {
            out.println(">" + index);
            out.println(ck.getKmerAsString());
            index++;
        }
    }
}
