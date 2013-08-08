package uk.ac.ox.well.indiana.analyses.reconstruction;

import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.fastq.FastqWriterFactory;
import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.*;
import uk.ac.ox.well.indiana.tools.Module;
import uk.ac.ox.well.indiana.utils.alignment.SmithWaterman;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

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
    //public PrintStream out;

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

    private SAMFileHeader createSAMFileHeader(Map<String, String> refMap) {
        SAMFileHeader samHeader = new SAMFileHeader();

        samHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        for (String geneName : refMap.keySet()) {
            SAMSequenceRecord samSeqRecord = new SAMSequenceRecord(geneName, refMap.get(geneName).length() + 1);
            samHeader.addSequence(samSeqRecord);
        }

        SAMReadGroupRecord contigReadGroup = new SAMReadGroupRecord("test");
        contigReadGroup.setSample("test");
        samHeader.addReadGroup(contigReadGroup);

        return samHeader;
    }

    private class SortedSW implements Comparable<SortedSW> {
        public SmithWaterman sw;
        public String contig;
        public String gene;

        public SortedSW(SmithWaterman sw, String contig, String gene) {
            this.sw = sw;
            this.contig = contig;
            this.gene = gene;
        }

        @Override
        public int compareTo(SortedSW o) {
            return Integer.valueOf(o.sw.getAlignmentScore()).compareTo(sw.getAlignmentScore());
        }

        public String toString() {
            return sw.getAlignmentScore() + " " + contig;
        }
    }

    private SortedSW alignRead(FastqRecord rec) {
        String fw = rec.getReadString();
        String rc = SequenceUtils.reverseComplement(fw);

        Set<SortedSW> alignments = new TreeSet<SortedSW>();

        for (String gene : targets.keySet()) {
            String target = targets.get(gene);

            SmithWaterman fsw = new SmithWaterman(target, fw);
            SmithWaterman rsw = new SmithWaterman(target, rc);

            SmithWaterman sw = (fsw.getAlignmentScore() > rsw.getAlignmentScore()) ? fsw : rsw;
            String contig    = (fsw.getAlignmentScore() > rsw.getAlignmentScore()) ? fw  : rc;

            SortedSW ssw = new SortedSW(sw, contig, gene);
            alignments.add(ssw);
        }

        return alignments.iterator().next();
    }

    private SAMFileWriter initializeSamWriter(SAMFileHeader sfh) {
//        SAMFileWriterFactory sfwf = new SAMFileWriterFactory();
//        sfwf.setCreateIndex(true);
//        SAMFileWriter sfw = sfwf.makeSAMOrBAMWriter(sfh, false, out);
//
//        return sfw;
        return null;
    }

    @Override
    public void execute() {
        loadGeneSequencesLookupTable();
        loadTargets();

        //SAMFileHeader sfh = createSAMFileHeader(targets);
        //SAMFileWriter sfw = initializeSamWriter(sfh);

        FastqWriterFactory fwf = new FastqWriterFactory();
        FastqWriter fqw1 = fwf.newWriter(fqo1);
        FastqWriter fqw2 = fwf.newWriter(fqo2);

        FastqReader f1 = new FastqReader(F1);
        FastqReader f2 = new FastqReader(F2);

        Iterator<FastqRecord> f1it = f1.iterator();
        Iterator<FastqRecord> f2it = f2.iterator();

        int numReadsSeen = 0;
        int numReadsFound = 0;
        while (f1it.hasNext() && f2it.hasNext()) { // && numReadsFound < 5) {
            numReadsSeen++;

            FastqRecord f1rec = f1it.next();
            FastqRecord f2rec = f2it.next();

            boolean f1recmaps = anyKmerMapsToROI(f1rec);
            boolean f2recmaps = anyKmerMapsToROI(f2rec);

            if (f1recmaps || f2recmaps) {
                numReadsFound++;

                fqw1.write(f1rec);
                fqw2.write(f2rec);

                /*
                SortedSW s1 = alignRead(f1rec);
                SortedSW s2 = alignRead(f2rec);

                if (!s1.gene.equalsIgnoreCase(s2.gene)) {
                    SAMReadGroupRecord rg = sfh.getReadGroup("test");
                    String readName = f1rec.getReadHeader().replaceAll("/1", "");

                    SAMRecord samRecord1 = new SAMRecord(sfh);
                    samRecord1.setReadBases(s1.contig.getBytes());
                    samRecord1.setReferenceName(s1.gene);
                    samRecord1.setAlignmentStart(s1.sw.getAlignmentStart() + 1);
                    samRecord1.setMappingQuality(99);
                    samRecord1.setCigar(s1.sw.getCigar());
                    samRecord1.setDuplicateReadFlag(false);
                    samRecord1.setReadPairedFlag(true);
                    samRecord1.setReadName(readName);
                    samRecord1.setFirstOfPairFlag(true);
                    samRecord1.setSecondOfPairFlag(false);
                    samRecord1.setReadNegativeStrandFlag(false);
                    samRecord1.setAttribute("RG", rg.getReadGroupId());
                    samRecord1.setMateUnmappedFlag(false);
                    samRecord1.setMateNegativeStrandFlag(false);
                    samRecord1.setMateAlignmentStart(s2.sw.getAlignmentStart() + 1);
                    samRecord1.setMateReferenceName(s2.gene);
                    sfw.addAlignment(samRecord1);

                    SAMRecord samRecord2 = new SAMRecord(sfh);
                    samRecord2.setReadBases(s2.contig.getBytes());
                    samRecord2.setReferenceName(s2.gene);
                    samRecord2.setAlignmentStart(s2.sw.getAlignmentStart() + 1);
                    samRecord2.setMappingQuality(99);
                    samRecord2.setCigar(s2.sw.getCigar());
                    samRecord2.setDuplicateReadFlag(false);
                    samRecord2.setReadPairedFlag(true);
                    samRecord2.setReadName(readName);
                    samRecord1.setFirstOfPairFlag(false);
                    samRecord1.setSecondOfPairFlag(true);
                    samRecord2.setReadNegativeStrandFlag(false);
                    samRecord2.setAttribute("RG", rg.getReadGroupId());
                    samRecord2.setMateUnmappedFlag(false);
                    samRecord2.setMateNegativeStrandFlag(false);
                    samRecord2.setMateAlignmentStart(s1.sw.getAlignmentStart() + 1);
                    samRecord2.setMateReferenceName(s1.gene);
                    sfw.addAlignment(samRecord2);
                }
                */
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
