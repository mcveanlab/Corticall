package uk.ac.ox.well.indiana.commands.caller.recover;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
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
    @Argument(fullName="graph", shortName="g", doc="Pedigree graph")
    public CortexGraph GRAPH;

    @Argument(fullName="dirty", shortName="d", doc="Dirty graph")
    public CortexGraph DIRTY;

    @Output
    public File out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(DIRTY.getSampleName(0));
        if (childColor < 0) {
            throw new IndianaException("Sample '" + DIRTY.getSampleName(0) + "' not found in pedigree graph");
        }

        log.info("Recovering removed kmers for color {} ({})", childColor, GRAPH.getSampleName(childColor));

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing records")
                .updateRecord(GRAPH.getNumRecords() / 10)
                .maxRecord(GRAPH.getNumRecords())
                .message("records processed")
                .make(log);

        int numRecordsRecovered = 0;

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(makeHeader(GRAPH.getHeader(), childColor));

        for (CortexRecord cr : GRAPH) {
            if (cr.getCoverage(childColor) > 0) {
                cgw.addRecord(cr);
            } else {
                int otherSamplesWithCoverage = 0;

                for (int c = 0; c < cr.getNumColors(); c++) {
                    if (c != childColor && cr.getCoverage(c) > 0) {
                        otherSamplesWithCoverage++;
                    }
                }

                if (otherSamplesWithCoverage > 0) {
                    CortexRecord dr = DIRTY.findRecord(cr.getCortexKmer());

                    if (dr != null && dr.getCoverage(0) > 0) {
                        long[] binaryKmer = new long[cr.getBinaryKmer().length];
                        int[] coverages = new int[cr.getCoverages().length];
                        byte[] edges = new byte[cr.getEdges().length];

                        System.arraycopy(cr.getBinaryKmer(), 0, binaryKmer, 0, binaryKmer.length);
                        System.arraycopy(cr.getCoverages(), 0, coverages, 0, coverages.length);
                        System.arraycopy(cr.getEdges(), 0, edges, 0, edges.length);

                        coverages[childColor] = dr.getCoverage(0);
                        edges[childColor] = dr.getEdges()[0];

                        int kmerSize = cr.getKmerSize();
                        int kmerBits = cr.getKmerBits();

                        CortexRecord nr = new CortexRecord(binaryKmer, coverages, edges, kmerSize, kmerBits);
                        cgw.addRecord(nr);

                        log.debug("old: {}", cr);
                        log.debug("new: {}", nr);
                        log.debug("---");

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
