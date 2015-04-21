package uk.ac.ox.well.indiana.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataTable;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class FindVariantsInContigs2 extends Module {
    @Argument(fullName="bed", shortName="bed", doc="Bed")
    public File BED;

    @Argument(fullName="bam", shortName="b", doc="BAM")
    public SAMFileReader BAM;

    @Output
    public File out;

    @Override
    public void execute() {
        IntervalTreeMap<Map<String, Object>> knownVariants = new IntervalTreeMap<Map<String, Object>>();

        TableReader tr = new TableReader(BED, new String[] { "chr", "start", "stop", "attrs" });

        Map<String, Boolean> variantWasSeen = new HashMap<String, Boolean>();
        Map<String, Map<String, Object>> variantInfo = new HashMap<String, Map<String, Object>>();

        log.info("Loading variation intervals...");
        for (Map<String, String> te : tr) {
            String chr = te.get("chr");
            int start = Integer.valueOf(te.get("start"));
            int stop = Integer.valueOf(te.get("stop"));

            Map<String, Object> attrs = new HashMap<String, Object>();
            for (String kv : te.get("attrs").split(";")) {
                String[] pair = kv.split("=");
                attrs.put(pair[0], pair[1]);
            }

            attrs.put("chr", chr);
            attrs.put("start", start);
            attrs.put("stop", stop);

            Interval variantInterval = new Interval(chr, start, stop);

            knownVariants.put(variantInterval, attrs);

            variantWasSeen.put((String) attrs.get("id"), false);
            variantInfo.put((String) attrs.get("id"), attrs);
        }
        log.info("  {} intervals", knownVariants.size());

        VariantContextWriterBuilder vcwb = new VariantContextWriterBuilder();
        vcwb.setOutputFile(out);
        vcwb.unsetOption(Options.INDEX_ON_THE_FLY);
        vcwb.setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
        VariantContextWriter vcw = vcwb.build();

        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        VCFHeader header = new VCFHeader(headerLines);

        vcw.writeHeader(header);

        log.info("Finding variants...");
        int numComplete = 0, numIncomplete = 0, numTotal = 0;
        int allelesMatch = 0, allelesMismatch = 0;
        for (SAMRecord contig : BAM) {
            Interval contigInterval = new Interval(contig.getReferenceName(), contig.getAlignmentStart(), contig.getAlignmentEnd());

            if (knownVariants.containsOverlapping(contigInterval)) {
                Collection<Map<String, Object>> bs = knownVariants.getOverlapping(contigInterval);

                log.debug("{} {}:{}-{}", contig.getReadName(), contig.getReferenceName(), contig.getAlignmentStart(), contig.getAlignmentEnd());
                for (Map<String, Object> b : bs) {
                    String id  = (String) b.get("id");
                    String ref = (String) b.get("ref");
                    String alt = (String) b.get("alt");

                    variantWasSeen.put(id, true);

                    if (alt != null) {
                        int bStart = (Integer) b.get("start");
                        int bStop = (Integer) b.get("stop");

                        int start = bStart - contig.getAlignmentStart();
                        int stop = bStop - contig.getAlignmentStart();
                        boolean isComplete = false;

                        if (start >= 0 && stop < contig.getReadLength()) {
                            isComplete = true;
                            numComplete++;
                        } else {
                            numIncomplete++;
                        }
                        numTotal++;

                        VariantContext newVC;

                        /*
                        b.remove("ref");
                        b.remove("alt");
                        b.remove("left");
                        b.remove("right");
                        */

                        if (contig.getReadNegativeStrandFlag()) {
                            newVC = (new VariantContextBuilder())
                                    .chr(contig.getReadName())
                                    .start(contig.getAlignmentEnd() - contig.getAlignmentStart() - stop)
                                    .stop(contig.getAlignmentEnd() - contig.getAlignmentStart() - start - 1)
                                    .alleles(SequenceUtils.reverseComplement(alt), SequenceUtils.reverseComplement(ref))
                                    .attributes(b)
                                    .attribute("isComplete", isComplete)
                                    .make();
                        } else {
                            newVC = (new VariantContextBuilder())
                                    .chr(contig.getReadName())
                                    .start(start + 1)
                                    .stop(stop)
                                    .alleles(alt, ref)
                                    .attributes(b)
                                    .attribute("isComplete", isComplete)
                                    .make();
                        }

                        String allele = "(incomplete)";
                        if (newVC.getAttributeAsBoolean("isComplete", false)) {
                            String seq = (contig.getReadNegativeStrandFlag()) ? SequenceUtils.reverseComplement(contig.getReadString()) : contig.getReadString();

                            allele = seq.substring(newVC.getStart(), newVC.getEnd() + 1);

                            if (!newVC.getReference().getBaseString().equals(allele)) {
                                log.debug("---");
                                log.debug("    - ref ({}): {}", contig.getReadNegativeStrandFlag() ? "rc" : "fw", newVC.getReference());
                                log.debug("    - sub ({}): {}", contig.getReadNegativeStrandFlag() ? "rc" : "fw", allele);
                                log.debug("---");

                                newVC = (new VariantContextBuilder(newVC))
                                        .filter("ALLELE_MISMATCH")
                                        .make();

                                allelesMismatch++;
                            } else {
                                allelesMatch++;
                            }
                        } else {
                            newVC = (new VariantContextBuilder(newVC))
                                    .filter("PARTIAL")
                                    .make();
                        }

                        log.debug("  {}-{}", bStart, bStop);
                        log.debug("      {}-{} {} {}", start, stop, isComplete, newVC);
                        if (contig.getReadNegativeStrandFlag()) {
                            log.debug("      ref (rc): {}", newVC.getReference());
                            log.debug("      sub (rc): {}", allele);
                        } else {
                            log.debug("      ref (fw): {}", newVC.getReference());
                            log.debug("      sub (fw): {}", allele);
                        }

                        vcw.add(newVC);
                    }
                }

                log.debug("");
            }
        }

        vcw.close();

        log.info("  complete={}/{}, incomplete={}/{}, total={}/{}", numComplete, numTotal, numIncomplete, numTotal, (numComplete + numIncomplete), numTotal);
        log.info("  match={} mismatch={}", allelesMatch, allelesMismatch);

        DataTable foundStats = new DataTable("foundStats", "found stats");
        foundStats.addColumns("type", "found", "missed");

        int numMissed = 0;
        for (String id : variantWasSeen.keySet()) {
            String type = (String) variantInfo.get(id).get("denovo");
            foundStats.set(type, "type", type);

            if (!variantWasSeen.get(id)) {
                log.info("  missed: {}", Joiner.on(", ").withKeyValueSeparator("=").join(variantInfo.get(id)));

                foundStats.increment(type, "missed");
            } else {
                foundStats.increment(type, "found");
            }
        }

        log.info("  missed={}", numMissed);

        log.info("\n{}", foundStats);
    }
}
