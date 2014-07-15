package uk.ac.ox.well.indiana.attic.analyses.LongContigs;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class FindSupernodesWithContigKmers extends Module {
    @Argument(fullName="supernodes", shortName="s", doc="Supernodes FASTA file")
    public FastaSequenceFile SUPERNODES;

    @Argument(fullName="contigs", shortName="c", doc="Contigs FASTA file")
    public ArrayList<FastaSequenceFile> CONTIGS;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Set<CortexKmer> kmers = new HashSet<CortexKmer>();

        for (FastaSequenceFile contigs : CONTIGS) {
            ReferenceSequence rseq;
            while ((rseq = contigs.nextSequence()) != null) {
                String seq = new String(rseq.getBases());

                for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                    CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                    kmers.add(kmer);
                }
            }
        }

        ReferenceSequence rseq;
        while ((rseq = SUPERNODES.nextSequence()) != null) {
            String seq = new String(rseq.getBases());

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                if (kmers.contains(kmer)) {
                    out.println(">" + rseq.getName() + "\n" + seq);
                    break;
                }
            }
        }
    }
}
