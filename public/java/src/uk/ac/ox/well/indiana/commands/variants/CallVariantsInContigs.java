package uk.ac.ox.well.indiana.commands.variants;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.File;
import java.util.*;

@Description(text="Call variants in contigs that have been aligned to a reference genome")
public class CallVariantsInContigs extends Module {
    @Argument(fullName="reference", shortName="r", doc="Samtools-indexed reference fasta file")
    public IndexedFastaSequenceFile REF;

    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public SAMFileReader BAM;

    @Argument(fullName="missingReadGroupSampleName", shortName="m", doc="Sample name to use for reads that lack a read group")
    public String MISSING_READ_GROUP_SAMPLE_NAME = "noReadGroup";

    @Output(fullName="out", shortName="o", doc="The output VCF file")
    public File out;

    private String getReferenceAllele(String chr, int pos) {
        return new String(REF.getSubsequenceAt(chr, pos, pos).getBases());
    }

    private Map<Integer, VariantContext> call(SAMRecord alignment) {
        StringBuilder target = new StringBuilder(new String(REF.getSubsequenceAt(alignment.getReferenceName(), alignment.getAlignmentStart(), alignment.getAlignmentEnd()).getBases()));
        StringBuilder query = new StringBuilder(alignment.getReadString());

        List<CigarElement> ces = alignment.getCigar().getCigarElements();
        CigarElement ceFirst = ces.get(0);
        CigarElement ceLast = ces.get(ces.size() - 1);

        if (ceLast.getOperator().equals(CigarOperator.S)) {
            query.delete(query.length() - ceLast.getLength(), query.length());
        }

        if (ceFirst.getOperator().equals(CigarOperator.S)) {
            query.delete(0, ceFirst.getLength());
        }

        String sampleName = MISSING_READ_GROUP_SAMPLE_NAME;
        if (alignment.getReadGroup() != null && alignment.getReadGroup().getSample() != null) {
            sampleName = alignment.getReadGroup().getSample();
        }

        Map<Integer, VariantContext> vcs = new TreeMap<Integer, VariantContext>();

        int pos = 0, tpos = 0, qpos = 0;
        for (CigarElement ce : alignment.getCigar().getCigarElements()) {
            if (ce.getOperator().equals(CigarOperator.M)) {
                for (int i = 0; i < ce.getLength(); i++) {
                    if (query.charAt(pos + i) != target.charAt(pos + i)) {
                        String targetAllele = String.valueOf(target.charAt(pos + i));
                        String queryAllele = String.valueOf(query.charAt(pos + i));

                        Genotype g = (new GenotypeBuilder(sampleName, Arrays.asList(Allele.create(queryAllele)))).make();

                        VariantContext vc = (new VariantContextBuilder())
                                .chr(alignment.getReferenceName())
                                .start(alignment.getAlignmentStart() + tpos + i)
                                .stop(alignment.getAlignmentStart() + tpos + i)
                                .alleles(targetAllele, queryAllele)
                                .attribute("contigName", alignment.getReadName())
                                .attribute("contigPos", qpos + i)
                                .attribute("mappingQuality", alignment.getMappingQuality())
                                .attribute("cigarString", alignment.getCigarString())
                                .attribute("alignmentLength", alignment.getAlignmentEnd() - alignment.getAlignmentStart())
                                .genotypes(g)
                                .make();

                        vcs.put(alignment.getAlignmentStart() + tpos + i, vc);
                    }
                }

                pos += ce.getLength();
                tpos += ce.getLength();
                qpos += ce.getLength();
            } else if (ce.getOperator().equals(CigarOperator.I)) {
                String targetAllele = pos == 0 ? getReferenceAllele(alignment.getReferenceName(), alignment.getAlignmentStart() - 1) : target.substring(pos - 1, pos);
                String queryAllele = pos == 0 ? getReferenceAllele(alignment.getReferenceName(), alignment.getAlignmentStart() - 1) + query.substring(pos, pos + ce.getLength()) : query.substring(pos - 1, pos + ce.getLength());

                Genotype g = (new GenotypeBuilder(sampleName, Arrays.asList(Allele.create(queryAllele)))).make();

                VariantContext vc = (new VariantContextBuilder())
                        .chr(alignment.getReferenceName())
                        .start(alignment.getAlignmentStart() + tpos - 1)
                        .stop(alignment.getAlignmentStart() + tpos - 1)
                        .alleles(targetAllele, queryAllele)
                        .attribute("contigName", alignment.getReadName())
                        .attribute("contigPos", qpos)
                        .attribute("mappingQuality", alignment.getMappingQuality())
                        .attribute("cigarString", alignment.getCigarString())
                        .attribute("alignmentLength", alignment.getAlignmentEnd() - alignment.getAlignmentStart())
                        .genotypes(g)
                        .make();

                vcs.put(alignment.getAlignmentStart() + tpos - 1, vc);

                target.insert(pos, StringUtil.repeatCharNTimes('-', ce.getLength()));
                pos += ce.getLength();
                qpos += ce.getLength();
            } else if (ce.getOperator().equals(CigarOperator.D)) {
                String targetAllele = pos == 0 ? getReferenceAllele(alignment.getReferenceName(), alignment.getAlignmentStart() - 1) + target.substring(pos + ce.getLength()) : target.substring(pos - 1, pos + ce.getLength());
                String queryAllele = pos == 0 ? getReferenceAllele(alignment.getReferenceName(), alignment.getAlignmentStart() - 1) : query.substring(pos - 1, pos);

                Genotype g = (new GenotypeBuilder(sampleName, Arrays.asList(Allele.create(queryAllele)))).make();

                VariantContext vc = (new VariantContextBuilder())
                        .chr(alignment.getReferenceName())
                        .start(alignment.getAlignmentStart() + tpos - 1)
                        .stop(alignment.getAlignmentStart() + tpos - 1 + ce.getLength())
                        .alleles(targetAllele, queryAllele)
                        .attribute("contigName", alignment.getReadName())
                        .attribute("contigPos", qpos)
                        .attribute("mappingQuality", alignment.getMappingQuality())
                        .attribute("cigarString", alignment.getCigarString())
                        .attribute("alignmentLength", alignment.getAlignmentEnd() - alignment.getAlignmentStart())
                        .genotypes(g)
                        .make();

                vcs.put(alignment.getAlignmentStart() + tpos - 1, vc);

                query.insert(pos, StringUtil.repeatCharNTimes('-', ce.getLength()));
                pos += ce.getLength();
                tpos += ce.getLength();
            }
        }

        if (vcs.size() > 0) {
            log.debug("samrecord: {}", alignment.getSAMString().trim());
            log.debug("t: {}", target);
            log.debug("q: {}", query);
            for (int i : vcs.keySet()) {
                log.debug("v: {}", vcs.get(i));
            }
            log.debug("-----");
        }

        return vcs;
    }

    @Override
    public void execute() {
        Map<String, Map<Integer, Set<VariantContext>>> variants = new TreeMap<String, Map<Integer, Set<VariantContext>>>();

        log.info("Calling variants...");
        Set<String> sampleNames = new HashSet<String>();
        int numRecords = 0;
        int numVariants = 0;
        int numContigsWithVariants = 0;

        for (SAMRecord sr : BAM) {
            if (numRecords % 5000 == 0) {
                log.info("  {} records", numRecords);
            }
            numRecords++;

            if (!sr.getReadUnmappedFlag()) {
                Map<Integer, VariantContext> vcs = call(sr);

                if (vcs.size() > 0) {
                    if (sr.getReadGroup() != null && sr.getReadGroup().getSample() != null) {
                        sampleNames.add(sr.getReadGroup().getSample());
                    } else {
                        sampleNames.add(MISSING_READ_GROUP_SAMPLE_NAME);
                    }
                }

                numVariants += vcs.size();
                numContigsWithVariants += (vcs.size() > 0) ? 1 : 0;

                for (int i : vcs.keySet()) {
                    String chrName = vcs.get(i).getChr();

                    if (!variants.containsKey(chrName)) {
                        variants.put(chrName, new TreeMap<Integer, Set<VariantContext>>());
                    }
                    if (!variants.get(chrName).containsKey(i)) {
                        variants.get(chrName).put(i, new HashSet<VariantContext>());
                    }

                    variants.get(chrName).get(i).add(vcs.get(i));
                }
            }
        }
        log.info("  found {} variants in {} contigs (some of these may be redundant)", numVariants, numContigsWithVariants);

        log.info("Writing VCF...");
        VariantContextWriterBuilder vcwb = new VariantContextWriterBuilder();
        vcwb.setOutputFile(out);
        vcwb.setReferenceDictionary(BAM.getFileHeader().getSequenceDictionary());
        VariantContextWriter vcw = vcwb.build();

        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        VCFStandardHeaderLines.addStandardFormatLines(headerLines, false, "GT");

        VCFHeader header = new VCFHeader(headerLines, sampleNames);
        header.setSequenceDictionary(BAM.getFileHeader().getSequenceDictionary());
        header.addMetaDataLine(new VCFInfoHeaderLine("contigName", 1, VCFHeaderLineType.String, "The name of the contig in which the variants were called"));
        header.addMetaDataLine(new VCFInfoHeaderLine("contigPos", 1, VCFHeaderLineType.Integer, "The 0-based position of the variant within the aligned portion of the contig"));
        header.addMetaDataLine(new VCFInfoHeaderLine("mappingQuality", 1, VCFHeaderLineType.Integer, "The mapping quality of the contig"));
        header.addMetaDataLine(new VCFInfoHeaderLine("cigarString", 1, VCFHeaderLineType.String, "The cigar string describing the alignment of the contig to the reference"));
        header.addMetaDataLine(new VCFInfoHeaderLine("alignmentLength", 1, VCFHeaderLineType.Integer, "The length of the aligned portion of the contig"));
        vcw.writeHeader(header);

        int numVariantsWritten = 0;
        for (SAMSequenceRecord ssr : BAM.getFileHeader().getSequenceDictionary().getSequences()) {
            String chrName = ssr.getSequenceName();
            for (int i : variants.get(chrName).keySet()) {
                for (VariantContext vc : variants.get(chrName).get(i)) {
                    vcw.add(vc);

                    numVariantsWritten++;
                }
            }
        }
        log.info("  wrote {} records", numVariantsWritten);

        vcw.close();
    }
}
