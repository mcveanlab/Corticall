package uk.ac.ox.well.cortexjdk.commands.discover.call;

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
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class CountNovelKmersInPartitions extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="roi", shortName="r", doc="ROIs")
    public CortexGraph ROI;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Set<CanonicalKmer> rois = new HashSet<>();
        for (CortexRecord cr : ROI) {
            rois.add(cr.getCanonicalKmer());
        }

        out.println(Joiner.on("\t").join("partitionName", "partitionLength", "novelKmers"));

        log.info("Processing partitions...");
        ReferenceSequence rseq;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            log.info("  {}", rseq.getName());

            Set<CanonicalKmer> cks = getUsedCanonicalKmers(rseq.getBaseString(), rois);

            out.println(Joiner.on("\t").join(rseq.getName().split(" ")[0], rseq.length(), cks.size()));
        }
    }

    private Set<CanonicalKmer> getUsedCanonicalKmers(String seq, Set<CanonicalKmer> rois) {
        Set<CanonicalKmer> used = new HashSet<>();

        for (int i = 0; i <= seq.length() - ROI.getKmerSize(); i++) {
            CanonicalKmer ck = new CanonicalKmer(seq.substring(i, i + ROI.getKmerSize()));

            if (rois.contains(ck)) {
                used.add(ck);
            }
        }

        return used;
    }
}
