package uk.ac.ox.well.indiana.attic.analyses.LongContigs;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.StringUtil;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;

import java.io.PrintStream;
import java.util.*;

public class KmerOrigin extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference")
    public IndexedFastaSequenceFile REFERENCE;

    @Argument(fullName="gff", shortName="g", doc="GFF")
    public GFF3 GFF;

    @Argument(fullName="sequence", shortName="s", doc="Query sequence")
    public HashMap<String, String> SEQUENCES;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CortexKmer, Set<String>> kmers = new HashMap<CortexKmer, Set<String>>();

        log.info("Processing reference...");

        ReferenceSequence rseq;
        while ((rseq = REFERENCE.nextSequence()) != null) {
            log.info("\t{}", rseq.getName());

            String seq = new String(rseq.getBases());

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                kmers.put(kmer, new HashSet<String>());
            }
        }

        log.info("Processing GFF...");

        for (GFF3Record gr : GFF) {
            if (gr.getType().equals("gene")) {
                rseq = REFERENCE.getSubsequenceAt(gr.getSeqid(), gr.getStart(), gr.getEnd());

                String seq = new String(rseq.getBases());

                for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                    CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                    kmers.get(kmer).add(gr.getAttribute("ID"));
                }
            }
        }

        for (String id : SEQUENCES.keySet()) {
            String seq = SEQUENCES.get(id);

            out.println(seq);

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                String padding = StringUtil.repeatCharNTimes(' ', i);

                Set<String> loci = kmers.get(kmer);
                for (String locus : loci) {
                    if (!locus.equals(id)) {
                        out.println(padding + locus);
                    }
                }
            }

            out.println();
        }
    }
}
