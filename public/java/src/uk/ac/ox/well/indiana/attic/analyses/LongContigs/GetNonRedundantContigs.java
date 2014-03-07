package uk.ac.ox.well.indiana.attic.analyses.LongContigs;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;

import java.io.PrintStream;
import java.util.*;

public class GetNonRedundantContigs extends Module {
    @Argument(fullName="fasta", shortName="f", doc="FASTA file")
    public FastaSequenceFile FASTA;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CortexKmer, Set<String>> kmerSharing = new HashMap<CortexKmer, Set<String>>();
        Map<String, String> seqs = new HashMap<String, String>();
        Map<String, Boolean> isContained = new HashMap<String, Boolean>();

        log.info("Loading sequences...");

        ReferenceSequence rseq;
        while ((rseq = FASTA.nextSequence()) != null) {
            String seq = new String(rseq.getBases());

            if (seq.length() > 0) {
                for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                    CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                    if (!kmerSharing.containsKey(kmer)) {
                        kmerSharing.put(kmer, new HashSet<String>());
                    }

                    kmerSharing.get(kmer).add(rseq.getName());
                }

                seqs.put(rseq.getName(), new String(rseq.getBases()));
                isContained.put(rseq.getName(), false);
            }
        }

        log.info("Identifying unique sequences...");

        int index = 0;
        for (CortexKmer kmer : kmerSharing.keySet()) {
            if (index % (kmerSharing.keySet().size() / 10) == 0) {
                log.info("  processed {}/{} (~{}%) kmers",
                         index,
                         kmerSharing.keySet().size(),
                         String.format("%.1f", 100.0*index/kmerSharing.keySet().size())
                );
            }
            index++;

            List<String> contigs = new ArrayList<String>();
            List<String> contigNames = new ArrayList<String>();

            for (String contigName : kmerSharing.get(kmer)) {
                contigs.add(seqs.get(contigName));
                contigNames.add(contigName);
            }

            for (int i = 0; i < contigs.size(); i++) {
                String c1name = contigNames.get(i);

                if (!isContained.get(c1name)) {
                    String c1 = contigs.get(i);
                    //String c1rc = SequenceUtils.reverseComplement(c1);

                    for (int j = 0; j < contigs.size(); j++) {
                        if (i != j) {
                            String c2 = contigs.get(j);

                            if (c1.contains(c2)) {
                                String c2name = contigNames.get(j);

                                isContained.put(c2name, true);
                            }
                        }
                    }
                }
            }
        }

        log.info("Writing unique sequences...");

        int numseqs = 0;
        for (String name : seqs.keySet()) {
            if (!isContained.get(name)) {
                String seq = seqs.get(name);

                out.println(">" + name);
                out.println(seq);

                numseqs++;
            }
        }

        log.info("  wrote {} non-redundant sequences", numseqs);
    }
}
