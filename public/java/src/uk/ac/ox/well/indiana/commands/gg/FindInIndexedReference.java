package uk.ac.ox.well.indiana.commands.gg;

import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.Set;

/**
 * Created by kiran on 28/12/2015.
 */
public class FindInIndexedReference extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference sequence")
    public File REFERENCE;

    @Argument(fullName="kmer", shortName="s", doc="Kmer sequence")
    public String KMER;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        KmerLookup kl = new KmerLookup(REFERENCE);

        String fw = KMER;
        String rc = SequenceUtils.reverseComplement(fw);

        Set<Interval> fwi = kl.findKmer(fw);
        Set<Interval> rci = kl.findKmer(rc);

        log.info("{} {}", fwi.size(), rci.size());

        for (Interval i : fwi) {
            log.info("fw {} {}", fw, i);
        }

        for (Interval i : rci) {
            log.info("rc {} {}", rc, i);
        }
    }
}
