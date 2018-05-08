package uk.ac.ox.well.cortexjdk.commands.discover.call;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;

public class ShowNovelKmers extends Module {
    @Argument(fullName = "contigs", shortName = "c", doc = "Contigs")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName = "roi", shortName = "r", doc = "ROI")
    public CortexGraph ROIS;

    @Argument(fullName = "graph", shortName = "g", doc = "Graph")
    public CortexGraph GRAPH;

    @Override
    public void execute() {
        ReferenceSequence rseq;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            String seq = rseq.getBaseString();

            log.info("{}", rseq.getName());

            for (int i = 0; i <= seq.length() - ROIS.getKmerSize(); i++) {
                String sk = seq.substring(i, i + ROIS.getKmerSize());
                CanonicalKmer ck = new CanonicalKmer(sk);

                log.info("{}/{} {} {} {}", i, seq.length() - ROIS.getKmerSize(), sk, ROIS.findRecord(ck) != null, GRAPH.findRecord(ck));
            }
        }
    }
}
