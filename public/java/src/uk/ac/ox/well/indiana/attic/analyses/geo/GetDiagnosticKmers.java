package uk.ac.ox.well.indiana.attic.analyses.geo;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.IntervalTree;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;

import java.io.PrintStream;
import java.util.*;

public class GetDiagnosticKmers extends Module {
    @Argument(fullName="reference", shortName="R", doc="Reference sequence")
    public IndexedFastaSequenceFile REFERENCE;

    @Argument(fullName="gff", shortName="gff", doc="GFF file")
    public GFF3 GFF;

    @Argument(fullName="sequenceName", shortName="sn", doc="Sequence names to include")
    public HashSet<String> SEQUENCES;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Set<CortexKmer> geneKmers = new HashSet<CortexKmer>();
        Map<String, IntervalTree<String>> geneRegions = new HashMap<String, IntervalTree<String>>();

        for (GFF3Record gr : GFF) {
            if (gr.getType().equals("gene")) {
                String id = gr.getAttribute("ID");

                if (SEQUENCES.contains(id)) {
                    log.info("Processing gene {}", id);

                    Collection<GFF3Record> exons = GFF3.getType("exon", GFF.getOverlapping(gr));

                    for (GFF3Record exon : exons) {
                        ReferenceSequence rseq = REFERENCE.getSubsequenceAt(exon.getSeqid(), exon.getStart(), exon.getEnd());
                        String seq = new String(rseq.getBases());

                        for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                            CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                            geneKmers.add(kmer);
                        }

                        if (!geneRegions.containsKey(gr.getSeqid())) {
                            geneRegions.put(gr.getSeqid(), new IntervalTree<String>());
                        }

                        geneRegions.get(gr.getSeqid()).put(exon.getStart(), exon.getEnd(), id);

                        //log.info("store name={}", gr.getSeqid());
                    }
                }
            }
        }

        log.info("Found {} gene kmers", geneKmers.size());

        Set<CortexKmer> kmerCopyNumber = new HashSet<CortexKmer>();

        ReferenceSequence rseq;
        while ((rseq = REFERENCE.nextSequence()) != null) {
            String[] names = rseq.getName().split("\\s+");
            String name = names[0];

            String seq = new String(rseq.getBases());

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                //log.info("recall name={}", name);
                if (geneRegions.containsKey(name) && geneRegions.get(name).minOverlapper(i, i+1) == null) {
                    CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                    kmerCopyNumber.add(kmer);
                }
            }
        }

        int numKmersFound = 0;
        for (CortexKmer kmer : geneKmers) {
            if (!kmerCopyNumber.contains(kmer)) {
                out.println(kmer);
                numKmersFound++;
            }
        }

        log.info("Found {} gene kmers unique to gene set (absent from rest of genome)", numKmersFound);
    }
}
