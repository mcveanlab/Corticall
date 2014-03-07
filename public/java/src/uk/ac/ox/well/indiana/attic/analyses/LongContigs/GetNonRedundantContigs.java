package uk.ac.ox.well.indiana.attic.analyses.LongContigs;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTreeMap;
import net.sf.samtools.*;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class GetNonRedundantContigs extends Module {
    @Argument(fullName="fasta", shortName="f", doc="FASTA file")
    public FastaSequenceFile FASTA;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, ReferenceSequence> seqs = new HashMap<String, ReferenceSequence>();
        Map<String, Boolean> isContained = new HashMap<String, Boolean>();

        log.info("Loading sequences...");
        ReferenceSequence rseq;
        while ((rseq = FASTA.nextSequence()) != null) {
            String seq = new String(rseq.getBases());

            if (seq.length() > 0) {
                seqs.put(seq, rseq);
                isContained.put(seq, false);
            }
        }

        log.info("Identifying unique sequences...");

        List<String> keys = new ArrayList<String>();
        keys.addAll(seqs.keySet());

        for (int i = 0; i < keys.size() - 1; i++) {
            if (i % (seqs.keySet().size() / 10) == 0) {
                log.info("  processed {}/{} sequences (~{}%)", i, seqs.keySet().size(), 100.0*i/seqs.keySet().size());
            }

            String seq1 = keys.get(i);
            for (int j = i + 1; j < keys.size(); j++) {
                String seq2 = keys.get(j);

                if (seq1.contains(seq2)) {
                    isContained.put(seq1, true);
                    break;
                }
            }
        }

        log.info("Writing unique sequences...");
        for (String seq1 : seqs.keySet()) {
            if (!isContained.get(seq1)) {
                ReferenceSequence rseq1 = seqs.get(seq1);

                out.println(">" + rseq1.getName());
                out.println(seq1);
            }
        }
    }
}
