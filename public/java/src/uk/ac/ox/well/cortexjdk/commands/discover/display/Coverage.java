package uk.ac.ox.well.cortexjdk.commands.discover.display;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;

import java.io.PrintStream;

public class Coverage extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="contigs", shortName="c", doc="Contigs")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="sample", shortName="s", doc="Sample")
    public String SAMPLE;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int c = GRAPH.getColorForSampleName(SAMPLE);

        out.println(Joiner.on("\t").join("contig", "kmer", "index", "coverage"));

        ReferenceSequence rseq;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            String seq = rseq.getBaseString();

            for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                String sk = seq.substring(i, i + GRAPH.getKmerSize());
                CanonicalKmer ck = new CanonicalKmer(sk);
                CortexRecord cr = GRAPH.findRecord(ck);
                int cov = cr.getCoverage(c);

                out.println(Joiner.on("\t").join(rseq.getName().split(" ")[0], sk, i, cov));
            }
        }
    }
}
