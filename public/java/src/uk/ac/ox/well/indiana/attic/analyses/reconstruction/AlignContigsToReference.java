package uk.ac.ox.well.indiana.attic.analyses.reconstruction;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.sw.SmithWaterman;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.util.*;

public class AlignContigsToReference extends Module {
    @Argument(fullName="contigTable", shortName="ct", doc="Contig table")
    public File CONTIG_TABLE;

    @Argument(fullName="reference", shortName="R", doc="Reference file")
    public FastaSequenceFile REFERENCE;

    @Output
    public File out;

    private class SAMComparator implements Comparator<SAMRecord> {
        @Override
        public int compare(SAMRecord o1, SAMRecord o2) {
            int s1 = o1.getIntegerAttribute("AL");
            int s2 = o2.getIntegerAttribute("AL");

            if (s1 == s2) { return 0; }
            return (s1 > s2) ? -1 : 1;
        }
    }

    private Map<String, String> loadReference() {
        Map<String, String> refMap = new TreeMap<String, String>();

        ReferenceSequence ref;
        while ((ref = REFERENCE.nextSequence()) != null) {
            String seq = new String(ref.getBases());

            String[] name = ref.getName().split("\\s+");

            refMap.put(name[0], seq);
        }

        return refMap;
    }

    private Set<String> getSampleList(TableReader tr) {
        Set<String> samples = new TreeSet<String>();
        for (Map<String, String> te : tr) {
            samples.add(te.get("sample"));
        }

        return samples;
    }

    private SAMFileHeader createSAMFileHeader(Map<String, String> refMap, Set<String> samples) {
        SAMFileHeader samHeader = new SAMFileHeader();

        samHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        for (String geneName : refMap.keySet()) {
            SAMSequenceRecord samSeqRecord = new SAMSequenceRecord(geneName, refMap.get(geneName).length() + 1);
            samHeader.addSequence(samSeqRecord);
        }

        for (String sampleName : samples) {
            SAMReadGroupRecord contigReadGroup = new SAMReadGroupRecord(sampleName + ".contigs");
            contigReadGroup.setSample(sampleName);
            samHeader.addReadGroup(contigReadGroup);

            SAMReadGroupRecord kmerReadGroup = new SAMReadGroupRecord(sampleName + ".unique_kmers");
            kmerReadGroup.setSample(sampleName);
            samHeader.addReadGroup(kmerReadGroup);
        }

        return samHeader;
    }

    @Override
    public void execute() {
        Map<String, String> refMap = loadReference();

        TableReader tr = new TableReader(CONTIG_TABLE);
        Set<String> samples = getSampleList(tr);

        SAMFileHeader sfh = createSAMFileHeader(refMap, samples);

        SAMFileWriterFactory sfwf = new SAMFileWriterFactory();
        sfwf.setCreateIndex(true);
        SAMFileWriter sfw = sfwf.makeSAMOrBAMWriter(sfh, false, out);

        int count = 0;
        for (Map<String, String> te : tr) {
            //if (count % (tr.size() / 10) == 0) {
                log.info("Processed {}/{} records", count, tr.size());
            //}
            count++;

            String sampleName = te.get("sample");
            List<String> kmerList = Arrays.asList(te.get("kmers").split(","));
            List<String> geneList = Arrays.asList(te.get("genes").split(","));
            Set<String> genes = new HashSet<String>(geneList);
            String fw = te.get("contig");
            String rc = SequenceUtils.reverseComplement(fw);

            for (int i = 0; i < kmerList.size(); i++) {
                String gene = geneList.get(i);
                String seq = refMap.get(gene);
                String kmerfw = kmerList.get(i);
                String kmerrc = SequenceUtils.reverseComplement(kmerfw);

                String kmer = null;

                if (seq.contains(kmerfw)) { kmer = kmerfw; }
                else if (seq.contains(kmerrc)) { kmer = kmerrc; }

                if (kmer != null) {
                    List<CigarElement> cigarElements = new ArrayList<CigarElement>();
                    cigarElements.add(new CigarElement(kmer.length(), CigarOperator.M));
                    Cigar cigar = new Cigar(cigarElements);

                    SAMReadGroupRecord rg = sfh.getReadGroup(sampleName + ".unique_kmers");

                    SAMRecord samRecord = new SAMRecord(sfh);
                    samRecord.setReadBases(kmer.getBytes());
                    samRecord.setAlignmentStart(seq.indexOf(kmer) + 1);
                    samRecord.setMappingQuality(0);
                    samRecord.setCigar(cigar);
                    samRecord.setDuplicateReadFlag(false);
                    samRecord.setReadPairedFlag(false);
                    samRecord.setReadName("km." + kmerfw.hashCode());
                    samRecord.setReferenceName(gene);
                    samRecord.setReadNegativeStrandFlag(false);
                    samRecord.setAttribute("RG", rg.getReadGroupId());
                    samRecord.setAttribute("SU", fw.hashCode());

                    sfw.addAlignment(samRecord);
                }

            }

            List<SAMRecord> reads = new ArrayList<SAMRecord>();

            for (String gene : genes) {
                String seq = refMap.get(gene);

                SmithWaterman fsw = new SmithWaterman(seq, fw);
                SmithWaterman rsw = new SmithWaterman(seq, rc);

                SmithWaterman sw = (fsw.getAlignmentScore() > rsw.getAlignmentScore()) ? fsw : rsw;
                String contig    = (fsw.getAlignmentScore() > rsw.getAlignmentScore()) ? fw  : rc;

                SAMReadGroupRecord rg = sfh.getReadGroup(sampleName + ".contigs");

                SAMRecord samRecord = new SAMRecord(sfh);
                samRecord.setReadBases(contig.getBytes());
                samRecord.setAlignmentStart(sw.getAlignmentStart() + 1);
                samRecord.setMappingQuality(99);
                samRecord.setCigar(sw.getCigar());
                samRecord.setDuplicateReadFlag(false);
                samRecord.setReadPairedFlag(false);
                samRecord.setReadName("cn." + fw.hashCode());
                samRecord.setReferenceName(gene);
                samRecord.setReadNegativeStrandFlag(false);
                samRecord.setAttribute("RG", rg.getReadGroupId());
                samRecord.setAttribute("SU", fw.hashCode());
                samRecord.setAttribute("AL", sw.getAlignmentScore());

                reads.add(samRecord);
            }

            Collections.sort(reads, new SAMComparator());

            for (int i = 0; i < reads.size(); i++) {
                if (i > 0) {
                    reads.get(i).setNotPrimaryAlignmentFlag(true);
                }

                if (reads.size() > 1) {
                    log.info("{} {} {} {} {}", sampleName, i, reads.get(i).getIntegerAttribute("AL"), reads.get(i).getReferenceName(), reads.get(i));
                }

                sfw.addAlignment(reads.get(i));
            }
        }

        sfw.close();
    }
}
