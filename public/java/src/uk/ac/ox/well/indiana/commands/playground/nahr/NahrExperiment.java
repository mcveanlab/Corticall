package uk.ac.ox.well.indiana.commands.playground.nahr;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.util.Arrays;
import java.util.List;

public class NahrExperiment extends Module {
    @Argument(fullName="reference", shortName="r", doc="Genome")
    public IndexedFastaSequenceFile REF;

    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Override
    public void execute() {
        Interval it1 = new Interval("Pf3D7_01_v3", 29510, 37126);   // PF3D7_0100100
        Interval it2 = new Interval("Pf3D7_02_v3", 916352, 923648); // PF3D7_0223500

        String var1 = REF.getSubsequenceAt(it1.getContig(), it1.getStart(), it1.getEnd()).getBaseString();
        String var2 = REF.getSubsequenceAt(it2.getContig(), it2.getStart(), it2.getEnd()).getBaseString();

        int childColor = GRAPH.getColorForSampleName("PG0063-C");
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(Arrays.asList("PG0051-C", "PG0052-C"));

        //DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs = CortexUtils.dfs(GRAPH, sk, childColor, parentColors, NahrStopper.class);

        for (int i = 0; i <= var1.length() - GRAPH.getKmerSize(); i++) {
            String sk = var1.substring(i, i + GRAPH.getKmerSize());
            CortexKmer ck = new CortexKmer(sk);

            CortexRecord cr = GRAPH.findRecord(ck);

            log.info("{}", cr);
        }
    }
}
