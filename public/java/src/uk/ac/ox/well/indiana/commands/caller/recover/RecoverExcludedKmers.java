package uk.ac.ox.well.indiana.commands.caller.recover;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraphWriter;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexHeader;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;

import java.io.File;

/**
 * Created by kiran on 20/06/2017.
 */
public class RecoverExcludedKmers extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="dirty", shortName="d", doc="Dirty graph")
    public CortexGraph DIRTY;

    @Argument(fullName="color", shortName="c", doc="Color")
    public Integer COLOR = 0;

    @Output
    public File out;

    @Override
    public void execute() {
        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing records")
                .updateRecord(GRAPH.getNumRecords() / 10)
                .message("records processed")
                .make(log);

        int numRecordsRecovered = 0;

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(makeHeader(GRAPH.getHeader(), COLOR));

        for (CortexRecord cr : GRAPH) {
            if (cr.getCoverage(COLOR) > 0) {
                cgw.addRecord(cr);
            } else {
                int otherSamplesWithCoverage = 0;

                for (int c = 0; c < cr.getNumColors(); c++) {
                    if (c != COLOR && cr.getCoverage(c) > 0) {
                        otherSamplesWithCoverage++;
                    }
                }

                //log.debug("{} {}", cr, otherSamplesWithCoverage);

                if (otherSamplesWithCoverage > 0) {
                    CortexRecord dr = DIRTY.findRecord(cr.getCortexKmer());

                    if (dr != null && dr.getCoverage(0) > 0) {
                        long[] binaryKmer = cr.getBinaryKmer();
                        int[] coverages = cr.getCoverages();
                        byte[] edges = cr.getEdges();
                        int kmerSize = cr.getKmerSize();
                        int kmerBits = cr.getKmerBits();

                        coverages[COLOR] = dr.getCoverage(0);
                        edges[COLOR] = dr.getEdges()[0];

                        CortexRecord nr = new CortexRecord(binaryKmer, coverages, edges, kmerSize, kmerBits);

                        cgw.addRecord(nr);

                        numRecordsRecovered++;
                    }
                }
            }

            pm.update();
        }

        log.info("Number of dirty records recovered: {}/{}", numRecordsRecovered, DIRTY.getNumRecords());
    }

    private CortexHeader makeHeader(CortexHeader fullHeader, int color) {
        CortexHeader newHeader = new CortexHeader();
        newHeader.setVersion(fullHeader.getVersion());
        newHeader.setKmerSize(fullHeader.getKmerSize());
        newHeader.setKmerBits(fullHeader.getKmerBits());
        newHeader.setNumColors(1);
        newHeader.addColor(fullHeader.getColor(color));

        return newHeader;
    }
}
