package uk.ac.ox.well.indiana.attic.analyses.geo;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IntervalTree;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.PrintStream;
import java.util.*;

public class GetDiagnosticKmers extends Module {
    @Argument(fullName="reference", shortName="R", doc="Reference sequence")
    public ArrayList<IndexedFastaSequenceFile> REFERENCES;

    @Argument(fullName="gff", shortName="gff", doc="GFF file")
    public ArrayList<GFF3> GFFS;

    @Argument(fullName="sequenceName", shortName="sn", doc="Sequence names to include", required=false)
    public HashSet<String> SEQUENCES;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Set<CortexKmer> geneKmers = new HashSet<CortexKmer>();
        Map<String, IntervalTree<String>> geneRegions = new HashMap<String, IntervalTree<String>>();

        log.info("Processing gene kmers...");
        for (int index = 0; index < REFERENCES.size(); index++) {
            IndexedFastaSequenceFile ref = REFERENCES.get(index);
            GFF3 gff = GFFS.get(index);

            for (GFF3Record gr : gff) {
                if (gr.getType().equals("gene")) {
                    String id = gr.getAttribute("ID");

                    if (SEQUENCES == null || SEQUENCES.isEmpty() || SEQUENCES.contains(id)) {
                        Collection<GFF3Record> exons = GFF3.getType("exon", gff.getOverlapping(gr));

                        for (GFF3Record exon : exons) {
                            ReferenceSequence rseq = ref.getSubsequenceAt(exon.getSeqid(), exon.getStart(), exon.getEnd());
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
        }

        log.info("\tfound {} gene kmers", geneKmers.size());

        log.info("Storing kmers outside of specified gene regions...");

        Set<CortexKmer> kmerCopyNumber = new HashSet<CortexKmer>();

        for (int index = 0; index < REFERENCES.size(); index++) {
            IndexedFastaSequenceFile ref = REFERENCES.get(index);

            ReferenceSequence rseq;
            while ((rseq = ref.nextSequence()) != null) {
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
        }

        log.info("\tfound {} genome kmers", kmerCopyNumber.size());

        log.info("Printing kmers...");

        TableWriter tw = new TableWriter(out);

        int numKmersFound = 0;
        for (CortexKmer kmer : geneKmers) {
            if (!kmerCopyNumber.contains(kmer)) {
                numKmersFound++;

                //out.println(">kmer." + numKmersFound);
                //out.println(kmer);

                Map<String, String> te = new LinkedHashMap<String, String>();
                te.put("gene", "multiple");
                te.put("kmer", kmer.getKmerAsString());

                tw.addEntry(te);
            }
        }

        log.info("\tfound {} gene kmers unique to gene set (absent from rest of genome)", numKmersFound);
    }
}
