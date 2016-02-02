package uk.ac.ox.well.indiana.commands.evaluate;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

public class RefBasedDenovoKmerStats extends Module {
    @Argument(fullName="denovoKmers", shortName="d", doc="Kmers spanning de novo events")
    public CortexGraph DENOVO_KMERS;

    @Argument(fullName="pacbioRef", shortName="p", doc="Kmers from PacBio ref")
    public CortexGraph PACBIO_KMERS;

    @Argument(fullName="cc", shortName="cc", doc="Child (clean)")
    public CortexGraph CHILD_CLEAN;

    @Argument(fullName="cd", shortName="cd", doc="Child (dirty)")
    public CortexGraph CHILD_DIRTY;

    @Argument(fullName="parents", shortName="pa", doc="Parents (dirty)")
    public CortexGraph PARENTS_DIRTY;

    @Override
    public void execute() {
        for (CortexRecord cr : DENOVO_KMERS) {
            CortexKmer ck = cr.getCortexKmer();

            CortexRecord crp  = PACBIO_KMERS.findRecord(ck);
            CortexRecord crcc = CHILD_CLEAN.findRecord(ck);
            CortexRecord crcd = CHILD_DIRTY.findRecord(ck);
            CortexRecord crpa = PARENTS_DIRTY.findRecord(ck);

            log.info("{}", ck);
        }
    }
}
