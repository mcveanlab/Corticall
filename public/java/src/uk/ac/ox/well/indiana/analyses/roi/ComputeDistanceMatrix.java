package uk.ac.ox.well.indiana.analyses.roi;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.tools.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataFrame;

import java.io.PrintStream;
import java.lang.ref.Reference;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class ComputeDistanceMatrix extends Module {
    @Argument(fullName="sequences", shortName="s", doc="Sequences in FASTA format")
    public FastaSequenceFile SEQUENCES;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Output
    public PrintStream out;

    public Set<String> getKmers(String seq, int kmerSize) {
        Set<String> kmers = new HashSet<String>();

        for (int i = 0; i <= seq.length() - kmerSize; i++) {
            String kmer = seq.substring(i, i + kmerSize);

            kmers.add(kmer);
        }

        return kmers;
    }

    public int kmerOverlap(Set<String> kmers1, Set<String> kmers2) {
        int overlap = 0;

        Set<String> allKmers = new HashSet<String>();
        allKmers.addAll(kmers1);
        allKmers.addAll(kmers2);

        for (String kmer : allKmers) {
            if (kmers1.contains(kmer) && kmers2.contains(kmer)) {
                overlap++;
            }
        }

        return overlap;
    }

    public int kmerTotal(Set<String> kmers1, Set<String> kmers2) {
        Set<String> allKmers = new HashSet<String>();
        allKmers.addAll(kmers1);
        allKmers.addAll(kmers2);

        return allKmers.size();
    }

    @Override
    public void execute() {
        Map<String, String> seqs = new HashMap<String, String>();

        ReferenceSequence rseq;
        while ((rseq = SEQUENCES.nextSequence()) != null) {
            seqs.put(rseq.getName(), new String(rseq.getBases()));
        }

        DataFrame<String, String, Float> dist = new DataFrame<String, String, Float>(-1.0f);

        for (String name1 : seqs.keySet()) {
            String seq1 = seqs.get(name1);
            Set<String> kmers1 = getKmers(seq1, KMER_SIZE);

            for (String name2 : seqs.keySet()) {
                if (dist.get(name1, name2) < 0) {
                    String seq2 = seqs.get(name2);
                    Set<String> kmers2 = getKmers(seq2, KMER_SIZE);

                    float overlap = kmerOverlap(kmers1, kmers2);
                    float total = kmerTotal(kmers1, kmers2);

                    dist.set(name1, name2, 1.0f - (overlap/total));
                    dist.set(name2, name1, 1.0f - (overlap/total));
                }
            }

        }

        out.println(dist);
    }
}
