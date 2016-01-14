package uk.ac.ox.well.indiana.commands.gg;

import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by kiran on 28/12/2015.
 */
public class FindInIndexedReference extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference sequence")
    public File REFERENCE;

    @Argument(fullName="kmer", shortName="s", doc="Kmer sequence", required=false)
    public String KMER;

    @Argument(fullName="novelKmers", shortName="n", doc="Novel kmers", required=false)
    public CortexGraph NOVEL_KMERS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        KmerLookup kl = new KmerLookup(REFERENCE);

        Set<String> kmersToCheck = new HashSet<String>();
        if (KMER != null) {
            kmersToCheck.add(KMER);
        }

        if (NOVEL_KMERS != null) {
            for (CortexRecord cr : NOVEL_KMERS) {
                kmersToCheck.add(cr.getKmerAsString());
            }
        }

        for (String fw : kmersToCheck) {
            String rc = SequenceUtils.reverseComplement(fw);

            Set<Interval> fwi = kl.findKmer(fw);
            Set<Interval> rci = kl.findKmer(rc);

            log.info("kmer: {} {} {}", fw, fwi.size(), rci.size());

            for (Interval i : fwi) {
                log.info("fw {} {}", fw, i);
            }

            for (Interval i : rci) {
                log.info("rc {} {}", rc, i);
            }
        }
    }
}
