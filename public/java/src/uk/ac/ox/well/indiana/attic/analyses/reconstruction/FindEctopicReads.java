package uk.ac.ox.well.indiana.attic.analyses.reconstruction;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class FindEctopicReads extends Module {
    @Argument(fullName="reference", shortName="R", doc="Reference sequence")
    public IndexedFastaSequenceFile FASTA;

    @Argument(fullName="fastq1", shortName="f1", doc="Fastq of pair first ends to process")
    public File F1;

    @Argument(fullName="fastq2", shortName="f2", doc="Fastq of pair second ends to process")
    public File F2;

    @Argument(fullName="gff", shortName="gff", doc="GFF3 file")
    public GFF3 GFF;

    @Argument(fullName="genes", shortName="g", doc="Genes to examine")
    public ArrayList<String> GENES;

    @Argument(fullName="window", shortName="w", doc="Window size")
    public Integer WINDOW = 500;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size for alignment")
    public Integer KMER_SIZE = 31;

    @Argument(fullName="targets", shortName="t", doc="Targets to align to")
    public FastaSequenceFile TARGETS;

    @Output(fullName="fastqOut1", shortName="o1", doc="Output for end 1")
    public File fqo1;

    @Output(fullName="fastqOut2", shortName="o2", doc="Output for end 2")
    public File fqo2;

    private class KmerInfo {
        public String gene;
        public String chr;
        public Integer pos;
    }

    private Map<CortexKmer, KmerInfo> geneSeqs;
    private Map<String, String> targets;

    private void loadTargets() {
        targets = new HashMap<String, String>();

        ReferenceSequence seq;
        while ((seq = TARGETS.nextSequence()) != null) {
            String sseq = new String(seq.getBases());
            String name = seq.getName();

            targets.put(name, sseq);
        }
    }

    private void loadGeneSequencesLookupTable() {
        geneSeqs = new HashMap<CortexKmer, KmerInfo>();

        for (String gene : GENES) {
            GFF3Record geneRecord = GFF.getRecord(gene);
            int start = geneRecord.getStart() - WINDOW < 0 ? 0 : geneRecord.getStart() - WINDOW;

            int seqLength = FASTA.getSequence(geneRecord.getSeqid()).length();
            int end = geneRecord.getEnd() + WINDOW > seqLength ? seqLength : geneRecord.getEnd() + WINDOW;

            String seq = new String(FASTA.getSubsequenceAt(geneRecord.getSeqid(), start, end).getBases());

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                KmerInfo ki = new KmerInfo();
                ki.gene = gene;
                ki.chr = geneRecord.getSeqid();
                ki.pos = start + i;

                geneSeqs.put(kmer, ki);
            }
        }
    }

    private boolean anyKmerMapsToROI(FastqRecord rec) {
        String seq = rec.getReadString();

        for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
            CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

            if (geneSeqs.containsKey(kmer)) {
                return true;
            }
        }

        return false;
    }

    @Override
    public void execute() {
        loadGeneSequencesLookupTable();
        loadTargets();

        FastqWriterFactory fwf = new FastqWriterFactory();
        FastqWriter fqw1 = fwf.newWriter(fqo1);
        FastqWriter fqw2 = fwf.newWriter(fqo2);

        FastqReader f1 = new FastqReader(F1);
        FastqReader f2 = new FastqReader(F2);

        Iterator<FastqRecord> f1it = f1.iterator();
        Iterator<FastqRecord> f2it = f2.iterator();

        int numReadsSeen = 0;
        int numReadsFound = 0;
        while (f1it.hasNext() && f2it.hasNext()) {
            numReadsSeen++;

            FastqRecord f1rec = f1it.next();
            FastqRecord f2rec = f2it.next();

            boolean f1recmaps = anyKmerMapsToROI(f1rec);
            boolean f2recmaps = anyKmerMapsToROI(f2rec);

            if (f1recmaps || f2recmaps) {
                numReadsFound++;

                fqw1.write(f1rec);
                fqw2.write(f2rec);
            }

            if (numReadsSeen % 1000000 == 0) {
                log.info("numReadsFound={}, numReadsSeen={}", numReadsFound, numReadsSeen);
            }
        }

        log.info("numReadsFound={}, numReadsSeen={}", numReadsFound, numReadsSeen);

        fqw1.close();
        fqw2.close();
    }
}
