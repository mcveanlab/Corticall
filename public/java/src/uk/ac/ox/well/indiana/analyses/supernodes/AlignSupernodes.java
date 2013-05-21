package uk.ac.ox.well.indiana.analyses.supernodes;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.*;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.alignment.SmithWaterman;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.util.*;

public class AlignSupernodes extends Tool {
    @Argument(fullName="referencePanel", shortName="rp", doc="Reference panel")
    public IndexedFastaSequenceFile REFERENCE;

    @Argument(fullName="supernodes", shortName="sn", doc="Related sequences file")
    public File SUPERNODES;

    @Output
    public File bamOut;

    private class SupernodeInfo {
        public String sample;
        public String supernode;
        public List<String> genes;
        //public List<SmithWaterman> alignments;
        //public List<Boolean> strand;
    }

    private ArrayList<SupernodeInfo> loadSupernodes() {
        ArrayList<SupernodeInfo> supernodes = new ArrayList<SupernodeInfo>();

        TableReader tr = new TableReader(SUPERNODES);
        for (Map<String, String> e : tr) {
            SupernodeInfo si = new SupernodeInfo();

            si.sample = e.get("sample");
            si.supernode = e.get("supernode");
            si.genes = new ArrayList<String>(Arrays.asList(e.get("genes").split(",")));

            supernodes.add(si);
        }

        return supernodes;
    }

    //          name    ref
    private Map<String, String> loadReferencePanel() {
        Map<String, String> panel = new HashMap<String, String>();

        ReferenceSequence ref;
        while ((ref = REFERENCE.nextSequence()) != null) {
            String seq = new String(ref.getBases());

            panel.put(ref.getName(), seq);
        }

        return panel;
    }

    @Override
    public void execute() {
        Map<String, String> panel = loadReferencePanel();

        SAMFileHeader samHeader = new SAMFileHeader();
        samHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        for (String gene : panel.keySet()) {
            if (!gene.contains(".")) {
                String seq = panel.get(gene);

                log.info("gene={}, length={}", gene, seq.length());

                samHeader.addSequence(new SAMSequenceRecord(gene, seq.length()));
            }
        }

        TableReader tr = new TableReader(SUPERNODES);
        Set<String> samples = new TreeSet<String>();
        for (Map<String, String> e : tr) {
            samples.add(e.get("sample"));
        }

        Map<String, String> rgIds = new HashMap<String, String>();
        for (String sample : samples) {
            SAMReadGroupRecord supernodeReadGroup = new SAMReadGroupRecord(sample + ".supernodes");
            supernodeReadGroup.setSample(sample);
            samHeader.addReadGroup(supernodeReadGroup);

            rgIds.put(sample, supernodeReadGroup.getReadGroupId());
        }

        SAMReadGroupRecord kmerReadGroup = new SAMReadGroupRecord("kmers");
        kmerReadGroup.setSample("kmers");
        samHeader.addReadGroup(kmerReadGroup);

        SAMFileWriterFactory sfwf = new SAMFileWriterFactory();
        sfwf.setCreateIndex(true);
        SAMFileWriter sfw = sfwf.makeSAMOrBAMWriter(samHeader, false, bamOut);

        for (String id : panel.keySet()) {
            if (id.contains(".")) {
                String kmer = panel.get(id);

                String[] pieces = id.split("\\.");
                String gene = pieces[0];

                String seq = panel.get(gene);
                int start = seq.indexOf(kmer);

                List<CigarElement> cigarElements = new ArrayList<CigarElement>();
                cigarElements.add(new CigarElement(kmer.length(), CigarOperator.M));
                Cigar cigar = new Cigar(cigarElements);

                SAMRecord kmerRecord = new SAMRecord(samHeader);
                kmerRecord.setReadBases(kmer.getBytes());
                kmerRecord.setAlignmentStart(start + 1);
                kmerRecord.setMappingQuality(0);
                kmerRecord.setCigar(cigar);
                kmerRecord.setDuplicateReadFlag(false);
                kmerRecord.setReadPairedFlag(false);
                kmerRecord.setReadName("km." + kmer.hashCode());
                kmerRecord.setReferenceName(gene);
                kmerRecord.setReadNegativeStrandFlag(false);
                kmerRecord.setAttribute("RG", kmerReadGroup.getReadGroupId());
                //kmerRecord.setAttribute("SU", supernode.hashCode());

                sfw.addAlignment(kmerRecord);
            }
        }

        tr = new TableReader(SUPERNODES);
        for (Map<String, String> e : tr) {
            String sample = e.get("sample");
            String supernode = e.get("supernode");
            String[] genes = e.get("genes").split(",");

            for (String gene : genes) {
                String seq = panel.get(gene);

                SmithWaterman fsw = new SmithWaterman(seq, supernode);
                SmithWaterman rsw = new SmithWaterman(seq, SequenceUtils.reverseComplement(supernode));

                SmithWaterman sw = (fsw.getAlignmentScore() > rsw.getAlignmentScore()) ? fsw : rsw;
                Boolean strand = (fsw.getAlignmentScore() > rsw.getAlignmentScore());

                String sn = (strand) ? supernode : SequenceUtils.reverseComplement(supernode);

                SAMRecord samRecord = new SAMRecord(samHeader);
                samRecord.setReadName("sn." + supernode.hashCode());
                samRecord.setReadBases(sn.getBytes());
                samRecord.setReferenceName(gene);
                samRecord.setAlignmentStart(sw.getAlignmentStart() + 1);
                samRecord.setCigar(sw.getCigar());
                samRecord.setMappingQuality(99);
                samRecord.setDuplicateReadFlag(false);
                samRecord.setReadPairedFlag(false);
                samRecord.setReadNegativeStrandFlag(false);
                samRecord.setAttribute("RG", rgIds.get(sample));
                samRecord.setAttribute("SU", supernode.hashCode());

                sfw.addAlignment(samRecord);

                if (genes.length > 1) {
                    //log.info("sn={} gene={} score={} start={} cigar={}", supernode.hashCode(), gene, sw.getAlignmentScore(), sw.getAlignmentStart(), sw.getCigar());
                    log.info("sn={} gene={} score={} start={}", supernode.hashCode(), gene, sw.getAlignmentScore(), sw.getAlignmentStart());
                }
            }
        }

        sfw.close();
    }
}
