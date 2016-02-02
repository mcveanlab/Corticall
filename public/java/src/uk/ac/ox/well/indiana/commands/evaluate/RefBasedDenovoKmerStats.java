package uk.ac.ox.well.indiana.commands.evaluate;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.PrintStream;
import java.util.LinkedHashMap;
import java.util.Map;

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

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TableWriter tw = new TableWriter(out);

        int records = 0;
        log.info("Processing ref-based novel kmers...");
        for (CortexRecord cr : DENOVO_KMERS) {
            if (records % (DENOVO_KMERS.getNumRecords()/10) == 0) {
                log.info("  {}/{} records", records, DENOVO_KMERS.getNumRecords());
            }

            CortexKmer ck = cr.getCortexKmer();

            boolean crp  = PACBIO_KMERS.findRecord(ck) != null;
            boolean crcc = CHILD_CLEAN.findRecord(ck) != null;
            boolean crcd = CHILD_DIRTY.findRecord(ck) != null;
            boolean crpa = PARENTS_DIRTY.findRecord(ck) != null;

            //log.info("{} inPacBio={} inClean={} inDirty={} inParents={}", ck, crp, crcc, crcd, crpa);

            Map<String, String> te = new LinkedHashMap<String, String>();
            te.put("kmer", ck.getKmerAsString());
            te.put("crp",  crp  ? "1" : "0");
            te.put("crcc", crcc ? "1" : "0");
            te.put("crcd", crcd ? "1" : "0");
            te.put("crpa", crpa ? "1" : "0");

            tw.addEntry(te);

            records++;
        }
    }
}
