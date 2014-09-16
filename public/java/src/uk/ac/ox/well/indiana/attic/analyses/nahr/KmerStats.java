package uk.ac.ox.well.indiana.attic.analyses.nahr;

import com.google.common.base.Joiner;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

import java.io.File;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class KmerStats extends Module {
    @Argument(fullName="graph", shortName="g", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="mask", shortName="m", doc="Cortex mask graph")
    public CortexGraph CORTEX_MASK;

    @Argument(fullName="covLimits", shortName="c", doc="Coverage limits file")
    public File COV_LIMITS;

    @Argument(fullName="sampleName", shortName="s", doc="Sample name")
    public String SAMPLE_NAME;

    @Argument(fullName="accession", shortName="a", doc="Accession")
    public String ACCESSION;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        log.info("Determining coverage threshold...");
        int covThreshold = 0;

        TableReader tr = new TableReader(COV_LIMITS);
        for (Map<String, String> te : tr) {
            if (te.get("sample").equals(SAMPLE_NAME)) {
                int p0Threshold = Integer.valueOf(te.get("p1"));
                int p1Threshold = Integer.valueOf(te.get("p2"));

                covThreshold = (p0Threshold > p1Threshold) ? p0Threshold : p1Threshold;
            }
        }

        log.info("  threshold = {}", covThreshold);

        log.info("Loading masked kmers...");
        Set<String> maskedKmers = new HashSet<String>();
        for (CortexRecord cr : CORTEX_MASK) {
            maskedKmers.add(cr.getKmerAsString());
        }

        log.info("  masked kmers: {}", maskedKmers.size());

        long rawKmers = CORTEX_GRAPH.getNumRecords();
        long excludedByMask = 0;
        long excludedByCoverage = 0;
        long remaining = 0;

        log.info("Examining raw kmers...");
        int numRecords = 0;
        for (CortexRecord cr : CORTEX_GRAPH) {
            if (numRecords % (CORTEX_GRAPH.getNumRecords() / 10) == 0) {
                log.info("  processed {}/{} records", numRecords, CORTEX_GRAPH.getNumRecords());
            }
            numRecords++;

            String kmer = cr.getKmerAsString();
            if (maskedKmers.contains(kmer)) {
                excludedByMask++;
            } else if (cr.getCoverage(0) <= covThreshold) {
                excludedByCoverage++;
            } else {
                remaining++;
            }
        }


        out.println(Joiner.on("\t").join("sampleName", "accession", "rawKmers", "excludedByMask", "excludedByCoverage", "remaining"));
        out.println(Joiner.on("\t").join(SAMPLE_NAME, ACCESSION, rawKmers, excludedByMask, excludedByCoverage, remaining));
    }
}
