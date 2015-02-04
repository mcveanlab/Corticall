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
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeMap;

public class FindVariantsInContigs extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VCF;

    @Argument(fullName="fastq", shortName="f", doc="Fastq flanks end 1")
    public File FASTQ;

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

            //log.info("  {}={}", name[0], seq.length());
        }
        log.info("  loaded {} sequences", contigSequences.size());

        log.info("Loading variants...");
        IntervalTreeMap<VariantContext> variants = new IntervalTreeMap<VariantContext>();
        for (VariantContext vc : VCF) {
            Interval interval = new Interval(vc.getChr(), vc.getStart(), vc.getStart());

            variants.put(interval, vc);
        }
        log.info("  loaded {} variants", variants.size());

        log.info("Loading fastq info...");
        Map<String, VariantContext> vinfo = new HashMap<String, VariantContext>();
        FastqReader fqr1 = new FastqReader(FASTQ);

        while (fqr1.hasNext()) {
            FastqRecord fqr = fqr1.next();

            String[] pieces = fqr.getReadHeader().split("[\\s:]+");
            String chr = pieces[1];
            int start = Integer.valueOf(pieces[2]);

            Interval interval = new Interval(chr, start, start);

            VariantContext vc = variants.get(interval);

            vinfo.put(pieces[0], vc);
        }
        log.info("  loaded {} records", vinfo.size());

        log.info("Loading flank alignments...");
        int numBothAligned = 0, numSameContig = 0, numUnambiguous = 0;
        for (SAMRecord sr : BAM) {
            VariantContext vc = vinfo.get(sr.getReadName());
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

            vinfo.put(sr.getReadName(), newvc);
        }
        log.info("  {} records with both flanks mapped", numBothAligned);
        log.info("  {} records with both flanks mapped and on the same contig", numSameContig);
        log.info("  {} records with both flanks mapped, on the same contig, with high MQ", numUnambiguous);

        log.info("Examining recovered contigs...");
        int numRecoveredVariants = 0;
        //Map<String, Integer> knownTypes = new HashMap<String, Integer>();
        //Map<String, Integer> recoveredTypes = new HashMap<String, Integer>();
        Map<String, Map<Integer, Integer>> knownTypes = new TreeMap<String, Map<Integer, Integer>>();
        Map<String, Map<Integer, Integer>> recoveredTypes = new TreeMap<String, Map<Integer, Integer>>();

        for (String id : vinfo.keySet()) {
            VariantContext newvc = vinfo.get(id);

            String type = newvc.getAttributeAsString("DENOVO", "unknown");
            int length = type.equals("DEL") ? newvc.getReference().length() - 1 : newvc.getAlternateAllele(0).length();

            if (type.equals("STR_CON")) {
                log.info("  str_con={}", newvc);
            }

            if (!knownTypes.containsKey(type)) {
                knownTypes.put(type, new TreeMap<Integer, Integer>());
            }

            if (!knownTypes.get(type).containsKey(length)) {
                knownTypes.get(type).put(length, 1);
            } else {
                knownTypes.get(type).put(length, knownTypes.get(type).get(length) + 1);
            }

            if (newvc.hasAttribute("left_flank_read") && newvc.hasAttribute("right_flank_read")) {
                SAMRecord end1 = (SAMRecord) newvc.getAttribute("left_flank_read");
                SAMRecord end2 = (SAMRecord) newvc.getAttribute("right_flank_read");

                if (end1.getReferenceName().equals(end2.getReferenceName()) &&
                    end1.getMappingQuality() > 0 && end2.getMappingQuality() > 0) {

                    String contig = contigSequences.get(end1.getReferenceName());

                    String leftFlankFw  = newvc.getAttributeAsString("flank_left", "N");
                    String altAlleleFw  = newvc.getAlternateAllele(0).getBaseString();
                    String rightFlankFw = newvc.getAttributeAsString("flank_right", "N");

                    String leftFlankRc  = SequenceUtils.reverseComplement(leftFlankFw);
                    String altAlleleRc  = SequenceUtils.reverseComplement(altAlleleFw);
                    String rightFlankRc = SequenceUtils.reverseComplement(rightFlankFw);

                    if ((contig.contains(leftFlankFw) && contig.contains(altAlleleFw) && contig.contains(rightFlankFw)) ||
                        (contig.contains(leftFlankRc) && contig.contains(altAlleleRc) && contig.contains(rightFlankRc))) {
                        numRecoveredVariants++;

                        if (!recoveredTypes.containsKey(type)) {
                            recoveredTypes.put(type, new TreeMap<Integer, Integer>());
                        }

                        if (!recoveredTypes.get(type).containsKey(length)) {
                            recoveredTypes.get(type).put(length, 1);
                        } else {
                            recoveredTypes.get(type).put(length, recoveredTypes.get(type).get(length) + 1);
                        }
                    }
                }
            }
        }
        log.info("  {} total recovered", numRecoveredVariants);
        //log.info("  knownTypes={}", Joiner.on(",").withKeyValueSeparator("=").join(knownTypes));
        //log.info("  recoveredTypes={}", Joiner.on(",").withKeyValueSeparator("=").join(recoveredTypes));

        TableWriter tw = new TableWriter(out);
        for (String type : knownTypes.keySet()) {
            for (int length : knownTypes.get(type).keySet()) {
                Map<String, String> te = new LinkedHashMap<String, String>();

                te.put("type", type);
                te.put("length", String.valueOf(length));
                te.put("known", String.valueOf(knownTypes.get(type).get(length)));
                te.put("recovered", String.valueOf(recoveredTypes.get(type).get(length)));

                tw.addEntry(te);
            }
        }
    }
}
