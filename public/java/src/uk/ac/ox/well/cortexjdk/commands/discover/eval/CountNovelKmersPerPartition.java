package uk.ac.ox.well.cortexjdk.commands.discover.eval;

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
import java.util.*;

/**
 * Created by kiran on 30/08/2017.
 */
public class CountNovelKmersPerPartition extends Module {
    @Argument(fullName="partitions", shortName="p", doc="Partitions")
    public FastaSequenceFile PARTITIONS;

    @Argument(fullName="rois", shortName="r", doc="ROIs")
    public CortexGraph ROIS;

    @Output
    public PrintStream out;

    @Output(fullName="rout", shortName="ro", doc="ROI info out")
    public PrintStream rout;

    @Override
    public void execute() {
        Map<CanonicalKmer, List<String>> rois = new HashMap<>();
        Map<CanonicalKmer, List<Integer>> rdist = new HashMap<>();

        for (CortexRecord cr : ROIS) {
            rois.put(cr.getCanonicalKmer(), new ArrayList<>());
            rdist.put(cr.getCanonicalKmer(), new ArrayList<>());
        }

        out.println(Joiner.on("\t").join("partition", "numNovels"));
        ReferenceSequence rseq;
        while ((rseq = PARTITIONS.nextSequence()) != null) {
            String rname = rseq.getName().split(" ")[0];
            String seq = rseq.getBaseString();

            int numNovels = 0;
            for (int i = 0; i <= seq.length() - ROIS.getKmerSize(); i++) {
                CanonicalKmer ck = new CanonicalKmer(seq.substring(i, i + ROIS.getKmerSize()));

                if (rois.containsKey(ck)) {
                    numNovels++;
                    rois.get(ck).add(rname);
                    rdist.get(ck).add(Math.min(i, seq.length() - ROIS.getKmerSize() - i));
                }
            }

            out.println(Joiner.on("\t").join(rname, numNovels));
        }

        rout.println(Joiner.on("\t").join("ck", "partition", "distanceToClosestEnd"));
        for (CanonicalKmer ck : rois.keySet()) {
            rout.println(Joiner.on("\t").join(
                    ck,
                    Joiner.on(",").join(rois.get(ck)),
                    Joiner.on(",").join(rdist.get(ck))
            ));
        }
    }
}
