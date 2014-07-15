package uk.ac.ox.well.indiana.attic.analyses.nahr.old;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;

import java.io.PrintStream;
import java.util.*;

public class FindEctopicRecombinationContigs extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="reference", shortName="r", doc="Reference")
    public IndexedFastaSequenceFile REF;

    @Argument(fullName="gff", shortName="g", doc="GFF")
    public GFF3 GFF;

    @Argument(fullName="gene1", shortName="g1", doc="Gene 1")
    public String GENE1;

    @Argument(fullName="gene2", shortName="g2", doc="Gene 2")
    public String GENE2;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Output
    public PrintStream out;

    private Map<CortexKmer, Integer> getGeneKmers(String gene) {
        Map<CortexKmer, Integer> geneKmers = new HashMap<CortexKmer, Integer>();

        GFF3Record gr = GFF.getRecord(gene);

        String seq = new String(REF.getSubsequenceAt(gr.getSeqid(), gr.getStart(), gr.getEnd()).getBases());
        for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
            CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

            if (!geneKmers.containsKey(kmer)) {
                geneKmers.put(kmer, 1);
            } else {
                geneKmers.put(kmer, geneKmers.get(kmer) + 1);
            }
        }

        return geneKmers;
    }

    @Override
    public void execute() {
        log.info("Processing gene kmers...");
        Map<CortexKmer, Integer> gene1Kmers = getGeneKmers(GENE1);
        Map<CortexKmer, Integer> gene2Kmers = getGeneKmers(GENE2);

        log.info("Processing reference kmers...");
        ReferenceSequence rseq;
        while ((rseq = REF.nextSequence()) != null) {
            String seq = new String(rseq.getBases());

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                if (gene1Kmers.containsKey(kmer)) {
                    gene1Kmers.put(kmer, gene1Kmers.get(kmer) + 1);
                }

                if (gene2Kmers.containsKey(kmer)) {
                    gene2Kmers.put(kmer, gene2Kmers.get(kmer) + 1);
                }
            }
        }

        log.info("Finding diagnostic kmers...");
        Set<CortexKmer> diagnostic1Kmers = new HashSet<CortexKmer>();
        for (CortexKmer kmer : gene1Kmers.keySet()) {
            if (gene1Kmers.get(kmer) == 2) {
                diagnostic1Kmers.add(kmer);
            }
        }

        Set<CortexKmer> diagnostic2Kmers = new HashSet<CortexKmer>();
        for (CortexKmer kmer : gene2Kmers.keySet()) {
            if (gene2Kmers.get(kmer) == 2) {
                diagnostic2Kmers.add(kmer);
            }
        }

        log.info("  found {}/{} diagnostic1 kmers and {}/{} diagnostic2 kmers", diagnostic1Kmers.size(), gene1Kmers.size(), diagnostic2Kmers.size(), gene2Kmers.size());

        log.info("Searching through contigs...");
        int contigsFound = 0;
        int contigsTotal = 0;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            String seq = new String(rseq.getBases());

            int g1 = 0, g2 = 0;

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                if (diagnostic1Kmers.contains(kmer)) { g1++; }
                if (diagnostic2Kmers.contains(kmer)) { g2++; }
            }

            if (g1 > 0 && g2 > 0) {
                out.println(">" + rseq.getName());
                out.println(seq);

                contigsFound++;
            }

            contigsTotal++;
        }

        log.info("  found {}/{} contigs", contigsFound, contigsTotal);
    }
}
