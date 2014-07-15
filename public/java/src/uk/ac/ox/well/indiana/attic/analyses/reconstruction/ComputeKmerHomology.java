package uk.ac.ox.well.indiana.attic.analyses.reconstruction;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;

import java.io.PrintStream;
import java.util.*;

public class ComputeKmerHomology extends Module {
    @Argument(fullName="reference", shortName="R", doc="Reference FASTA")
    public IndexedFastaSequenceFile REFERENCE;

    @Argument(fullName="gff", shortName="gff", doc="GFF file")
    public GFF3 GFF;

    @Argument(fullName="genes", shortName="g", doc="Genes to process")
    public ArrayList<String> GENES;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Output
    public PrintStream out;

    private class KmerInfo {
        public KmerInfo(String gene) {
            this.gene = gene;
            this.isPresent = false;
        }

        public String gene;
        public Boolean isPresent;
    }

    @Override
    public void execute() {
        //Map<CortexKmer, KmerInfo> kmers = new HashMap<CortexKmer, KmerInfo>();
        Map<String, Set<CortexKmer>> kmers = new HashMap<String, Set<CortexKmer>>();

        for (String gene : GENES) {
            GFF3Record gr = GFF.getRecord(gene);

            if (gr != null) {
                String chr = gr.getSeqid();
                int start = gr.getStart();
                int end = gr.getEnd();

                String seq = new String(REFERENCE.getSubsequenceAt(chr, start, end).getBases());
                for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                    CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                    //kmers.put(kmer, new KmerInfo(gene));

                    if (!kmers.containsKey(gene)) {
                        kmers.put(gene, new HashSet<CortexKmer>());
                    }

                    kmers.get(gene).add(kmer);
                }
            }
        }

        Map<String, Map<String, Integer>> kmerSharingAmongGenes = new TreeMap<String, Map<String, Integer>>();

        for (String gene1 : kmers.keySet()) {
            Set<CortexKmer> geneKmers = kmers.get(gene1);

            for (String gene2 : kmers.keySet()) {
                int count = 0;
                for (CortexKmer kmer : geneKmers) {
                    if (kmers.get(gene2).contains(kmer)) {
                        count++;
                    }
                }

                if (!kmerSharingAmongGenes.containsKey(gene1)) {
                    kmerSharingAmongGenes.put(gene1, new TreeMap<String, Integer>());
                }

                kmerSharingAmongGenes.get(gene1).put(gene2, count);
            }
        }

        List<String> names = new ArrayList<String>();
        for (String gene1 : kmerSharingAmongGenes.keySet()) {
            names.add(gene1);
        }

        out.println("\t" + Joiner.on("\t").join(names));

        for (String gene1 : kmerSharingAmongGenes.keySet()) {
            List<Integer> counts = new ArrayList<Integer>();
            for (String gene2 : kmerSharingAmongGenes.get(gene1).keySet()) {
                int count = kmerSharingAmongGenes.get(gene1).get(gene2);

                counts.add(count);
            }

            out.println(gene1 + "\t" + Joiner.on("\t").join(counts));
        }
    }
}
