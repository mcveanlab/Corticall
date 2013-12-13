package uk.ac.ox.well.indiana.analyses.geo;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.tools.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataFrame;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class ComputeSequenceRelatedness extends Module {
    @Argument(fullName="sequences", shortName="s", doc="Sequences")
    public FastaSequenceFile FASTA;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 11;

    @Argument(fullName="seqNames", shortName="sn", doc="Seq names", required=false)
    public HashSet<String> SEQ_NAMES;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, Set<String>> kmerToSeqMap = new HashMap<String, Set<String>>();

        ReferenceSequence rseq;
        while ((rseq = FASTA.nextSequence()) != null) {
            String[] names = rseq.getName().split("\\s+");
            String name = names[0];

            if (SEQ_NAMES == null || SEQ_NAMES.size() == 0 || SEQ_NAMES.contains(name)) {
                String bases = new String(rseq.getBases());

                for (int i = 0; i <= bases.length() - KMER_SIZE; i++) {
                    String kmer = bases.substring(i, i + KMER_SIZE);

                    if (!kmerToSeqMap.containsKey(kmer)) {
                        kmerToSeqMap.put(kmer, new HashSet<String>());
                    }

                    kmerToSeqMap.get(kmer).add(name);
                }
            }
        }

        DataFrame<String, String, Float> rmat = new DataFrame<String, String, Float>(0.0f);

        Map<String, Float> totalKmers = new HashMap<String, Float>();

        for (String kmer : kmerToSeqMap.keySet()) {
            for (String name1 : kmerToSeqMap.get(kmer)) {
                for (String name2 : kmerToSeqMap.get(kmer)) {
                    rmat.set(name1, name2, rmat.get(name1, name2) + 1.0f);
                }

                if (!totalKmers.containsKey(name1)) {
                    totalKmers.put(name1, 0.0f);
                }
                totalKmers.put(name1, totalKmers.get(name1) + 1.0f);
            }
        }

        for (String name1 : rmat.getRowNames()) {
            for (String name2 : rmat.getColNames()) {
                float intersection = rmat.get(name1, name2);
                float total = totalKmers.get(name1) + totalKmers.get(name2) - intersection;
                //rmat.set(name1, name2, 1.0f - (intersection/total));
                rmat.set(name1, name2, 1.0f - (intersection/total));
            }
        }

        out.println(rmat);
    }
}
