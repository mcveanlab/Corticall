package uk.ac.ox.well.cortexjdk.commands.discover.call;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertexFactory;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalUtils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class TrimPartitions extends Module {
    @Argument(fullName="partitions", shortName="p", doc="Partitions")
    public FastaSequenceFile FASTA;

    @Argument(fullName="rois", shortName="r", doc="ROIs")
    public CortexGraph ROI;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Set<CanonicalKmer> rois = new HashSet<>();
        for (CortexRecord rr : ROI) {
            rois.add(rr.getCanonicalKmer());
        }

        ReferenceSequence rseq;
        while ((rseq = FASTA.nextSequence()) != null) {
            String seq = rseq.getBaseString();

            List<CortexVertex> w = new ArrayList<>();
            int start = Integer.MAX_VALUE;
            int stop = 0;
            for (int i = 0; i <= seq.length() - ROI.getKmerSize(); i++) {
                String sk = seq.substring(i, i + ROI.getKmerSize());
                CanonicalKmer ck = new CanonicalKmer(sk);

                if (rois.contains(ck)) {
                    if (i < start) { start = i; }
                    if (i > stop) { stop = i; }
                }

                w.add(new CortexVertexFactory().bases(sk).make());
            }

            start = (start - 100 >= 0) ? start - 100 : 0;
            stop = (stop + 100 < w.size() - 1) ? stop + 100 : w.size() - 1;

            out.println(">" + rseq.getName());
            out.println(TraversalUtils.toContig(w.subList(start, stop)));
        }
    }
}
