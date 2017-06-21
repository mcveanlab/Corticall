package uk.ac.ox.well.indiana.commands.caller.call;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

/**
 * Created by kiran on 21/06/2017.
 */
public class AnnotateNahrCandidates extends Module {
    @Argument(fullName="sequences", shortName="s", doc="Sequences")
    public FastaSequenceFile SEQUENCES;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="refs", shortName="R", doc="References")
    public HashMap<String, KmerLookup> LOOKUPS;

    @Override
    public void execute() {
        List<String> backgrounds = new ArrayList<>(LOOKUPS.keySet());

        ReferenceSequence rseq;
        while ((rseq = SEQUENCES.nextSequence()) != null) {
            String seq = rseq.getBaseString();

            log.info("{}", rseq.getName());

            for (int i = 0; i <= seq.length() - ROI.getKmerSize(); i++) {
                String sk = seq.substring(i, i + ROI.getKmerSize());

                List<String> intervals = new ArrayList<>();

                for (String background : backgrounds) {
                    Set<Interval> its = LOOKUPS.get(background).findKmer(sk);
                    if (its.size() == 1) {
                        Interval it = its.iterator().next();
                        intervals.add(background + ":" + it.getContig() + ":" + it.getStart() + "-" + it.getEnd() + ":" + (it.isPositiveStrand() ? "+" : "-"));
                    } else {
                        intervals.add(background + ":none");
                    }
                }

                log.info("  {} {} {}", i, sk, Joiner.on("\t").join(intervals));
            }
        }
    }
}
