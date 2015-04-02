package uk.ac.ox.well.indiana.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.collections.comparators.ReverseComparator;
import org.apache.commons.lang.StringUtils;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class FindVariantsInContigs extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VCF;

    @Argument(fullName="bam", shortName="b", doc="BAM")
    public SAMFileReader BAM;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        log.info("Loading contigs...");
        Map<String, String> contigSequences = new HashMap<String, String>();
        ReferenceSequence rseq;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            String[] name = rseq.getName().split("\\s+");
            String seq = new String(rseq.getBases());

            contigSequences.put(name[0], seq);
        }
        log.info("  loaded {} sequences", contigSequences.size());

        log.info("Loading variants...");
        Map<String, VariantContext> variants = new HashMap<String, VariantContext>();
        for (VariantContext vc : VCF) {
            variants.put(vc.getAttributeAsString("SIMID", "unknown"), vc);
        }
        log.info("  loaded {} variants", variants.size());

        log.info("Loading flank alignments...");
        int numBothAligned = 0, numSameContig = 0, numUnambiguous = 0;
        for (SAMRecord sr : BAM) {
            VariantContext vc = variants.get(sr.getReadName());
            VariantContext newvc;

            if (sr.getFirstOfPairFlag()) {
                newvc = (new VariantContextBuilder(vc))
                        .attribute("left_flank_read", sr)
                        .make();
            } else {
                newvc = (new VariantContextBuilder(vc))
                        .attribute("right_flank_read", sr)
                        .make();
            }

            if (newvc.hasAttribute("left_flank_read") && newvc.hasAttribute("right_flank_read")) {
                numBothAligned++;

                SAMRecord end1 = (SAMRecord) newvc.getAttribute("left_flank_read");
                SAMRecord end2 = (SAMRecord) newvc.getAttribute("right_flank_read");

                if (end1.getReferenceName().equals(end2.getReferenceName())) {
                    numSameContig++;

                    if (end1.getMappingQuality() > 0 && end2.getMappingQuality() > 0) {
                        numUnambiguous++;
                    }
                }
            }

            variants.put(sr.getReadName(), newvc);
        }
        log.info("  {} records with both flanks mapped", numBothAligned);
        log.info("  {} records with both flanks mapped and on the same contig", numSameContig);
        log.info("  {} records with both flanks mapped, on the same contig, with high MQ", numUnambiguous);

        log.info("Examining recovered contigs...");
        Map<String, Integer> knownTypes = new HashMap<String, Integer>();
        TableWriter tw = new TableWriter(out);

        Map<String, Set<VariantContext>> gcVariants = new TreeMap<String, Set<VariantContext>>();
        Set<String> idsToRemove = new HashSet<String>();
        for (String id : variants.keySet()) {
            VariantContext newvc = variants.get(id);

            String gcid = newvc.getAttributeAsString("GCINDEX", "NA");
            String nahr = newvc.getAttributeAsString("NAHR", "NA");
            String denovo = newvc.getAttributeAsString("DENOVO", "NA");

            if (!gcid.equals("NA")) {
                String totalId = newvc.getChr() + ".gcid" + gcid;

                if (!gcVariants.containsKey(totalId)) {
                    gcVariants.put(totalId, new HashSet<VariantContext>());
                }

                gcVariants.get(totalId).add(newvc);

                idsToRemove.add(id);
            } else if (!nahr.equals("NA")) {
                idsToRemove.add(id);
            } else if (denovo.equals("RECOMB")) {
                idsToRemove.add(id);
            }
        }

        for (String totalId : gcVariants.keySet()) {
            //log.info("  id={}", totalId);

            for (VariantContext newvc : gcVariants.get(totalId)) {
                //log.info("    vc={}", newvc);
            }
        }

        for (String id : idsToRemove) {
            variants.remove(id);
        }

        for (String id : variants.keySet()) {
            VariantContext newvc = variants.get(id);

            String denovo = newvc.getAttributeAsString("DENOVO", "NA");
            String nahr = newvc.getAttributeAsString("NAHR", "NA");
            String event = nahr.equals("NA") ? denovo : "NAHR";
            String type = newvc.getType().name();
            String gcid = newvc.getAttributeAsString("GCINDEX", "NA");
            String repUnit = newvc.getAttributeAsString("REPUNIT", "NA");
            int repsBefore = newvc.getAttributeAsInt("REPS_BEFORE", 0);
            int repsAfter = newvc.getAttributeAsInt("REPS_AFTER", 0);
            boolean variantFound = false;
            String matchedSeq = "NA";
            int matchedPos = 0;
            int variantStart = 0;
            int variantEnd = 0;
            int pos1 = -1;
            int pos2 = -1;
            int length = 0;
            String flank1 = "NA";
            String flank2 = "NA";
            String refAllele = "NA";
            String altAllele = "NA";

            if (newvc.isSimpleInsertion()) {
                length = newvc.getAlternateAllele(0).length() - newvc.getReference().length();
            } else if (newvc.isSimpleDeletion()) {
                length = newvc.getReference().length() - newvc.getAlternateAllele(0).length();
            } else if (newvc.isMNP()) {
                length = newvc.getReference().length();
            }

            SAMRecord end1 = newvc.hasAttribute("left_flank_read") ? (SAMRecord) newvc.getAttribute("left_flank_read") : null;
            SAMRecord end2 = newvc.hasAttribute("right_flank_read") ? (SAMRecord) newvc.getAttribute("right_flank_read") : null;

            if (end1 != null && end2 != null &&
                end1.getReferenceName().equals(end2.getReferenceName()) &&
                end1.getMappingQuality() > 0 && end2.getMappingQuality() > 0) {

                String contig = contigSequences.get(end1.getReferenceName());

                String leftFlankFw  = newvc.getAttributeAsString("flank_left", "N");
                String refAlleleFw  = newvc.getReference().getBaseString();
                String altAlleleFw  = newvc.isVariant() ? newvc.getAlternateAllele(0).getBaseString() : newvc.getReference().getBaseString();
                String rightFlankFw = newvc.getAttributeAsString("flank_right", "N");
                String patternFw    = ".*(" + leftFlankFw + ")" + altAlleleFw + "(" + rightFlankFw + ").*";

                String leftFlankRc  = SequenceUtils.reverseComplement(leftFlankFw);
                String refAlleleRc  = SequenceUtils.reverseComplement(refAlleleFw);
                String altAlleleRc  = SequenceUtils.reverseComplement(altAlleleFw);
                String rightFlankRc = SequenceUtils.reverseComplement(rightFlankFw);
                String patternRc    = ".*(" + rightFlankRc + ")" + altAlleleRc + "(" + leftFlankRc + ").*";

                if (contig.contains(leftFlankFw) && contig.contains(rightFlankFw)) {
                    Pattern pFw = Pattern.compile(patternFw);
                    Matcher mFw = pFw.matcher(contig);

                    if (mFw.matches()) {
                        matchedSeq = mFw.group(0);
                        matchedPos = mFw.start(1);
                        pos1 = mFw.start(1);
                        pos2 = mFw.start(2);

                        variantFound = true;

                        flank1 = mFw.group(1);
                        flank2 = mFw.group(2);
                        refAllele = refAlleleFw;
                        altAllele = altAlleleFw;

                        variantStart = mFw.end(1) + 1;
                        variantEnd = mFw.start(2);
                    }
                } else if (contig.contains(leftFlankRc) && contig.contains(rightFlankRc)) {
                    Pattern pRc = Pattern.compile(patternRc);
                    Matcher mRc = pRc.matcher(contig);

                    if (mRc.matches()) {
                        matchedSeq = mRc.group(0);
                        matchedPos = mRc.start(1);
                        pos1 = mRc.start(1);
                        pos2 = mRc.start(2);

                        variantFound = true;

                        flank1 = mRc.group(1);
                        flank2 = mRc.group(2);
                        refAllele = refAlleleRc;
                        altAllele = altAlleleRc;

                        variantStart = mRc.end(1) + 1;
                        variantEnd = mRc.start(2);
                    }
                }
            }

            Map<String, String> te = new LinkedHashMap<String, String>();
            te.put("id", id);
            te.put("event", event);
            te.put("type", type);
            te.put("length", String.valueOf(length));
            te.put("variantFound", variantFound ? "TRUE" : "FALSE");
            te.put("matchedPos", String.valueOf(matchedPos));
            te.put("variantStart", String.valueOf(variantStart));
            te.put("variantEnd", String.valueOf(variantEnd));
            te.put("chr", newvc.getChr());
            te.put("start", String.valueOf(newvc.getStart()));
            te.put("gcindex", gcid);
            te.put("nahr", nahr);
            te.put("repUnit", repUnit);
            te.put("repsBefore", String.valueOf(repsBefore));
            te.put("repsAfter", String.valueOf(repsAfter));
            te.put("flank1", flank1);
            te.put("flank2", flank2);
            te.put("mq1", end1 == null ? "NA" : String.valueOf(end1.getMappingQuality()));
            te.put("mq2", end2 == null ? "NA" : String.valueOf(end2.getMappingQuality()));
            te.put("contig1", end1 == null ? "NA" : end1.getReferenceName());
            te.put("contig2", end2 == null ? "NA" : end2.getReferenceName());
            te.put("pos1", pos1 < 0 ? "NA" : String.valueOf(pos1));
            te.put("pos2", pos2 < 0 ? "NA" : String.valueOf(pos2));
            te.put("matchedSeq", matchedSeq);
            te.put("refAllele", refAllele);
            te.put("altAllele", altAllele);

            if (!event.equals("NA")) {
                tw.addEntry(te);

                if (!knownTypes.containsKey(event)) {
                    knownTypes.put(event, 1);
                } else {
                    knownTypes.put(event, knownTypes.get(event) + 1);
                }
            }
        }

        log.info("  {}", Joiner.on(", ").withKeyValueSeparator("=").join(knownTypes));
    }
}
