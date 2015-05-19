package uk.ac.ox.well.indiana.commands.gg;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import scala.ref.Reference;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class GetPathableContigs extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference")
    public FastaSequenceFile REFERENCE;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Output
    public PrintStream out;

    public void increment(Map<CortexKmer, Integer> kmerCopyNumber, CortexKmer ck) {
        if (!kmerCopyNumber.containsKey(ck)) {
            kmerCopyNumber.put(ck, 1);
        } else {
            kmerCopyNumber.put(ck, kmerCopyNumber.get(ck) + 1);
        }
    }

    @Override
    public void execute() {
        Map<CortexKmer, Integer> kmerCopyNumber = new HashMap<CortexKmer, Integer>();

        log.info("Processing contigs...");
        List<ReferenceSequence> chrs = new ArrayList<ReferenceSequence>();

        ReferenceSequence rseq;
        while ((rseq = REFERENCE.nextSequence()) != null) {
            log.info("  {}", rseq.getName());

            String seq = new String(rseq.getBases());

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                CortexKmer ck = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                increment(kmerCopyNumber, ck);
            }

            chrs.add(rseq);
        }

        log.info("Emitting contigs...");
        int numContigs = 0;
        for (ReferenceSequence cseq : chrs) {
            log.info("  {}", cseq.getName());

            String chr = new String(cseq.getBases());
            String[] name = cseq.getName().split("\\s+");

            for (int pos = 1; pos <= chr.length() - KMER_SIZE - 1; pos++) {
                CortexKmer ck = new CortexKmer(chr.substring(pos, pos + KMER_SIZE));

                if (kmerCopyNumber.get(ck) > 1) {
                    int start = pos;
                    CortexKmer pk = ck;
                    do {
                        start--;
                        pk = new CortexKmer(chr.substring(start, start + KMER_SIZE));
                    } while (start > 0 && kmerCopyNumber.get(pk) != 1);

                    int end = pos;
                    CortexKmer nk = ck;
                    while (end - start < 200 && end < chr.length() - KMER_SIZE - 1) {
                        do {
                            end++;
                            nk = new CortexKmer(chr.substring(end, end + KMER_SIZE));
                        } while (end < chr.length() - KMER_SIZE - 1 && kmerCopyNumber.get(nk) != 1);
                    }

                    String contig = chr.substring(start, end);
                    out.println(">" + name[0] + ".pos" + start + "-" + end + ".contig" + numContigs);
                    out.println(contig);

                    numContigs++;

                    pos = end;
                }
            }
        }

        log.info("  emitted {} contigs", numContigs);
    }
}
