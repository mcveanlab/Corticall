package uk.ac.ox.well.cortexjdk.commands.utils;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Description;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexCollection;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraphWriter;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

@Description(text = "Remove from primary graph kmers that appear in secondary graphs")
public class Remove extends Module {
    @Argument(fullName="graph", shortName="g", doc="Primary graph")
    public CortexGraph PGRAPH;

    @Argument(fullName="secondary", shortName="s", doc="Secondary graph")
    public ArrayList<CortexGraph> SGRAPH;

    @Output
    public File out;

    @Override
    public void execute() {
        List<CortexGraph> graphs = new ArrayList<>();
        graphs.add(PGRAPH);
        graphs.addAll(SGRAPH);
        CortexCollection cc = new CortexCollection(graphs);

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(PGRAPH.getHeader());

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing records...")
                .message("processed")
                .updateRecord(1000000)
                .make(log);

        long numKept = 0, numRemoved = 0;

        for (CortexRecord cr : cc) {
            boolean found = false;
            for (int c = PGRAPH.getNumColors(); c < cr.getNumColors(); c++) {
                if (cr.getCoverage(c) > 0) {
                    found = true;
                    break;
                }
            }

            if (!found) {
                int[] cov = new int[PGRAPH.getNumColors()];
                byte[] edges = new byte[PGRAPH.getNumColors()];

                for (int c = 0; c < PGRAPH.getNumColors(); c++) {
                    cov[c] = cr.getCoverage(c);
                    edges[c] = cr.getEdges()[c];
                }

                CortexRecord ncr = new CortexRecord(
                        cr.getBinaryKmer(),
                        cov,
                        edges,
                        cr.getKmerSize(),
                        cr.getKmerBits()
                );

                cgw.addRecord(ncr);

                numKept++;
            } else {
                numRemoved++;
            }

            pm.update("processed (kept: " + numKept + ", removed: " + numRemoved + ")");
        }

        log.info("  finished (kept: {}, removed: {})", numKept, numRemoved);

        cgw.close();
    }
}
