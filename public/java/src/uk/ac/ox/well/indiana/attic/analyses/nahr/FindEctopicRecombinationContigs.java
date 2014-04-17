package uk.ac.ox.well.indiana.attic.analyses.nahr;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
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

    @Argument(fullName="genes", shortName="genes", doc="Genes")
    public ArrayList<String> GENES;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CortexKmer, Integer> geneKmers = new HashMap<CortexKmer, Integer>();

        log.info("Processing gene kmers...");
        for (String gene : GENES) {
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
        }

        log.info("Processing reference kmers...");
        ReferenceSequence rseq;
        while ((rseq = REF.nextSequence()) != null) {
            String seq = new String(rseq.getBases());

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                if (geneKmers.containsKey(kmer)) {
                    geneKmers.put(kmer, geneKmers.get(kmer) + 1);
                }
            }
        }

        log.info("Finding diagnostic kmers...");
        Set<CortexKmer> diagnosticKmers = new HashSet<CortexKmer>();
        for (CortexKmer kmer : geneKmers.keySet()) {
            if (geneKmers.get(kmer) == 2) {
                diagnosticKmers.add(kmer);
            }
        }
        log.info("  found {}/{} diagnostic kmers", diagnosticKmers.size(), geneKmers.size());

        log.info("Searching through contigs...");
        int contigsFound = 0;
        int contigsTotal = 0;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            String seq = new String(rseq.getBases());
            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                if (diagnosticKmers.contains(kmer)) {
                    out.println(">" + rseq.getName());
                    out.println(seq);

                    contigsFound++;
                    break;
                }
            }

            contigsTotal++;
        }

        log.info("  found {}/{} contigs", contigsFound, contigsTotal);
    }
}
