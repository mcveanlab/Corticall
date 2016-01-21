package uk.ac.ox.well.indiana.commands.gg;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.PrintStream;
import java.util.HashSet;
import java.util.Set;

public class DiscardOrphans extends Module {
    @Argument(fullName="clean", shortName="c", doc="Graph")
    public CortexGraph CLEAN;

    @Argument(fullName="novelKmers", shortName="n", doc="Novel kmers")
    public CortexGraph NOVEL_KMERS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Set<CortexKmer> novelKmersToDiscard = new HashSet<CortexKmer>();

        log.info("Examining novel kmers for orphaned stretches...");

        for (CortexRecord cr : NOVEL_KMERS) {
            CortexKmer ck = cr.getCortexKmer();
            String sk = ck.getKmerAsString();

            DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs = CortexUtils.dfs(CLEAN, null, sk, 0, null, ChildTraversalStopper.class);

            Set<CortexKmer> novelKmers = new HashSet<CortexKmer>();
            boolean isOrphaned = true;

            for (AnnotatedVertex av : dfs.vertexSet()) {
                CortexKmer ak = new CortexKmer(av.getKmer());
                CortexRecord ar = CLEAN.findRecord(ak);

                if (ar != null) {
                    if (ar.getCoverage(0) > 0 && ar.getCoverage(1) == 0 && ar.getCoverage(2) == 0) {
                        novelKmers.add(ak);
                    } else if (ar.getCoverage(1) > 0 || ar.getCoverage(2) > 0) {
                        isOrphaned = false;
                        break;
                    }
                }
            }

            if (isOrphaned) {
                novelKmersToDiscard.addAll(novelKmers);
            }
        }

        log.info("  {} kmers to discard", novelKmersToDiscard.size());

        int index = 0;
        for (CortexKmer ck : novelKmersToDiscard) {
            out.println(">" + index);
            out.println(ck);
            index++;
        }
    }
}
