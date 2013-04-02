package uk.ac.ox.well.indiana.analyses.kmerSharing;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.*;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.alignment.SmithWaterman;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.io.utils.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class ExportSupernodesToBAM extends Tool {
    @Argument(fullName="reference", shortName="R", doc="Reference sequence")
    public IndexedFastaSequenceFile REFERENCE;

    @Argument(fullName="gff", shortName="gff", doc="GFF file")
    public GFF3 GFF;

    @Argument(fullName="relatedSequences", shortName="rel", doc="Related sequences file")
    public File RELATED_SEQUENCES;

    @Argument(fullName="sample", shortName="sn", doc="Sample name to place in resultant BAM file")
    public String SAMPLE_NAME;

    @Output
    public PrintStream out;

    @Output(fullName="bamOut", shortName="bo", doc="Output file for BAM")
    public File bamOut;

    @Output(fullName="fastaOut", shortName="fo", doc="Output file for FASTA")
    public PrintStream fastaOut;

    @Output(fullName="gffOut", shortName="go", doc="Output file for GFF")
    public PrintStream gffOut;

    @Override
    public int execute() {
        TreeMap<String, TreeMap<String, TreeSet<String>>> relatedSequences = loadRelatedSequences(RELATED_SEQUENCES);

        gffOut.println("##gff-version\t3");
        gffOut.println("##feature-ontology\tso.obo");
        gffOut.println("##attribute-ontology\tgff3_attributes.obo");

        SAMFileHeader samHeader = new SAMFileHeader();
        for (String gene : relatedSequences.keySet()) {
            GFF3Record record = GFF.getRecord(gene);
            String seq = new String(REFERENCE.getSubsequenceAt(record.getSeqid(), record.getStart(), record.getEnd()).getBases());

            SAMSequenceRecord samSeqRecord = new SAMSequenceRecord(gene, seq.length());

            samHeader.addSequence(samSeqRecord);

            fastaOut.println(">" + gene);
            fastaOut.println(seq);

            gffOut.println("##sequence-region\t" + gene + "\t1\t" + seq.length());
        }

        for (String gene : relatedSequences.keySet()) {
            GFF3Record geneRecord = GFF.getRecord(gene);

            Collection<GFF3Record> records = GFF.getOverlapping(geneRecord);

            for (GFF3Record record : records) {
                int start = record.getStart() - geneRecord.getStart() + 1;
                int end = start + (record.getEnd() - record.getStart());

                GFF3Record newRecord = new GFF3Record(record);

                newRecord.setSeqid(gene);
                newRecord.setStart(start);
                newRecord.setEnd(end);

                gffOut.println(newRecord);
            }
        }

        samHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        SAMReadGroupRecord supernodeReadGroup = new SAMReadGroupRecord(SAMPLE_NAME + ".supernodes");
        supernodeReadGroup.setSample(SAMPLE_NAME);
        samHeader.addReadGroup(supernodeReadGroup);

        SAMReadGroupRecord kmerReadGroup = new SAMReadGroupRecord(SAMPLE_NAME + ".unique_kmers");
        kmerReadGroup.setSample(SAMPLE_NAME);
        samHeader.addReadGroup(kmerReadGroup);

        SAMFileWriterFactory sfwf = new SAMFileWriterFactory();
        sfwf.setCreateIndex(true);
        SAMFileWriter sfw = sfwf.makeSAMOrBAMWriter(samHeader, false, bamOut);

        HashSet<Integer> seenSupernodes = new HashSet<Integer>();

        for (String gene : relatedSequences.keySet()) {
            log.info("Processing supernodes for gene '{}'", gene);

            GFF3Record record = GFF.getRecord(gene);
            String seq = new String(REFERENCE.getSubsequenceAt(record.getSeqid(), record.getStart(), record.getEnd()).getBases());

            TreeMap<String, TreeSet<String>> supernodesAndKmers = relatedSequences.get(gene);

            for (String fwsupernode : supernodesAndKmers.keySet()) {
                String rcsupernode = SequenceUtils.getReverseComplement(fwsupernode);

                SmithWaterman fsw = new SmithWaterman(seq, fwsupernode);
                SmithWaterman rsw = new SmithWaterman(seq, rcsupernode);

                SmithWaterman sw = fsw.getAlignmentScore() > rsw.getAlignmentScore() ? fsw : rsw;
                String supernode = fsw.getAlignmentScore() > rsw.getAlignmentScore() ? fwsupernode : rcsupernode;

                if (!seenSupernodes.contains(supernode.hashCode())) {
                    SAMRecord samRecord = new SAMRecord(samHeader);
                    samRecord.setReadBases(supernode.getBytes());
                    samRecord.setAlignmentStart(sw.getAlignmentStart() + 1);
                    samRecord.setMappingQuality(99);
                    samRecord.setCigar(sw.getCigar());
                    samRecord.setDuplicateReadFlag(false);
                    samRecord.setReadPairedFlag(false);
                    samRecord.setReadName("sn." + supernode.hashCode());
                    samRecord.setReferenceName(gene);
                    samRecord.setReadNegativeStrandFlag(false);
                    samRecord.setAttribute("RG", supernodeReadGroup.getReadGroupId());
                    samRecord.setAttribute("SU", supernode.hashCode());

                    byte[] baseQualities = new byte[supernode.length()];
                    for (int i = 0; i < baseQualities.length; i++) {
                        baseQualities[i] = 40;
                    }

                    TreeSet<String> kmers = supernodesAndKmers.get(fwsupernode);
                    for (String fwkmer : kmers) {
                        String rckmer = SequenceUtils.getReverseComplement(fwkmer);

                        String kmer = seq.contains(fwkmer) ? fwkmer : rckmer;
                        int start = seq.contains(fwkmer) ? seq.indexOf(fwkmer) : seq.indexOf(rckmer);

                        if (supernode.contains(kmer)) {
                            int kmerIndex = supernode.indexOf(kmer);
                            for (int i = kmerIndex; i < kmerIndex + kmer.length(); i++) {
                                baseQualities[i] = 0;
                            }
                        }

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
                        kmerRecord.setReadName("sn." + supernode.hashCode() + ".km." + kmer.hashCode());
                        kmerRecord.setReferenceName(gene);
                        kmerRecord.setReadNegativeStrandFlag(false);
                        kmerRecord.setAttribute("RG", kmerReadGroup.getReadGroupId());
                        kmerRecord.setAttribute("SU", supernode.hashCode());

                        sfw.addAlignment(kmerRecord);
                    }

                    samRecord.setBaseQualities(baseQualities);
                    sfw.addAlignment(samRecord);

                    seenSupernodes.add(supernode.hashCode());
                }
            }
        }

        sfw.close();

        return 0;
    }

    private TreeMap<String, TreeMap<String, TreeSet<String>>> loadRelatedSequences(File relatedSequenceFile) {
        //      gene            supernode       kmer
        TreeMap<String, TreeMap<String, TreeSet<String>>> relatedSequences = new TreeMap<String, TreeMap<String, TreeSet<String>>>();

        TableReader reader = new TableReader(relatedSequenceFile);
        for (HashMap<String, String> entry : reader) {
            String gene = entry.get("genes");
            String kmer = entry.get("kmer");
            String supernode = entry.get("superNode");

            if (!relatedSequences.containsKey(gene)) {
                relatedSequences.put(gene, new TreeMap<String, TreeSet<String>>());
            }

            if (!relatedSequences.get(gene).containsKey(supernode)) {
                relatedSequences.get(gene).put(supernode, new TreeSet<String>());
            }

            relatedSequences.get(gene).get(supernode).add(kmer);
        }

        return relatedSequences;
    }
}
