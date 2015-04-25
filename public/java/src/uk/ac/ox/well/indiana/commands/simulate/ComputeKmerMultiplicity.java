package uk.ac.ox.well.indiana.commands.simulate;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;

import java.io.File;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Set;

public class ComputeKmerMultiplicity extends Module {
    @Argument(fullName="kmerList", shortName="l", doc="List of kmers")
    public File KMER_LIST;

    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Set<CortexKmer> kmers = new HashSet<CortexKmer>();
        LineReader lr = new LineReader(KMER_LIST);
        while (lr.hasNext()) {
            String l = lr.getNextRecord();

            kmers.add(new CortexKmer(l));
        }

        for (CortexRecord cr : GRAPH) {
            CortexKmer ck = cr.getCortexKmer();

            if (kmers.contains(ck)) {
                out.println(ck.getKmerAsString() + "\t" + cr.getCoverage(0));
            }
        }
    }
}
