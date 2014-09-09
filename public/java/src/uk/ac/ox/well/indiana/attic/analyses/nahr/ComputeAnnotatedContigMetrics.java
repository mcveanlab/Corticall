package uk.ac.ox.well.indiana.attic.analyses.nahr;

import com.google.common.base.Joiner;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class ComputeAnnotatedContigMetrics extends Module {
    @Argument(fullName="annotatedContigs", shortName="ac", doc="Annotated contigs")
    public ArrayList<File> ANNS;

    @Argument(fullName="bamsCanonical", shortName="bc", doc="BAMs (canonical)")
    public ArrayList<SAMFileReader> BAMSC;

    @Argument(fullName="bams0", shortName="b0", doc="BAMs (for parent 0)")
    public ArrayList<SAMFileReader> BAMS0;

    @Argument(fullName="bams1", shortName="b0", doc="BAMs (for parent 1)")
    public ArrayList<SAMFileReader> BAMS1;

    @Argument(fullName="deltas", shortName="d", doc="Delta")
    public ArrayList<File> DELTAS;

    @Argument(fullName="gff", shortName="g", doc="GFF file")
    public GFF3 GFF;

    //@Argument(fullName="recombs", shortName="rc", doc="Recombs info")
    //public File RECOMBS;

    @Output
    public PrintStream out;

    private class ContigInfo {
        public Interval canonicalLocus;
        public Interval ref0Locus;
        public Interval ref1Locus;
        public Interval ref0ToCanonicalBlock;
        public Interval ref1ToCanonicalBlock;
        public Interval ref0ToCanonicalExact;
        public Interval ref1ToCanonicalExact;
        public boolean sameChromosome = false;
        public boolean sameEffectiveLocus = false;
        public boolean alignedToCanonical = false;
        public boolean alignedToRef0 = false;
        public boolean alignedToRef1 = false;
        public boolean clippedInCanonical = false;
        public boolean clippedInRef0 = false;
        public boolean clippedInRef1 = false;
        public boolean perfectAlignmentInCanonical = false;
        public boolean perfectAlignmentInRef0 = false;
        public boolean perfectAlignmentInRef1 = false;
        public int numAlignmentsInCanonical = 0;
        public int numAlignmentsInRef0 = 0;
        public int numAlignmentsInRef1 = 0;
    }

    private int longestRun(String ann, char entry) {
        Set<Integer> runs = new HashSet<Integer>();

        int l = 0;
        for (int i = 0; i < ann.length(); i++) {
            if (ann.charAt(i) == entry) {
                l++;
            } else {
                runs.add(l);

                l = 0;
            }
        }

        int longestRunLength = 0;
        for (int runLength : runs) {
            if (runLength > longestRunLength) {
                longestRunLength = runLength;
            }
        }

        return longestRunLength;
    }

    private int numSwitches(String ann) {
        String ann2 = ann.replaceAll("B", "").replaceAll("\\.", "");
        int numSwitches = 0;

        if (ann2.length() > 0) {
            char c = ann2.charAt(0);

            for (int i = 0; i < ann2.length(); i++) {
                if (ann2.charAt(i) != c) {
                    numSwitches++;
                    c = ann2.charAt(i);
                }
            }
        }

        return numSwitches;
    }

    private IntervalTreeMap<Interval> loadDeltas(boolean reverse) {
        IntervalTreeMap<Interval> imap = new IntervalTreeMap<Interval>();

        for (File delta : DELTAS) {
            LineReader lr = new LineReader(delta);

            String currentHeader = null;
            String refContig = null;
            String altContig = null;

            String line;
            while ((line = lr.getNextRecord()) != null) {
                if (line.startsWith(">")) {
                    currentHeader = line;

                    line = line.replace(">", "");
                    String[] fields = line.split("\\s+");
                    refContig = fields[0];
                    altContig = fields[1];
                }
                else {
                    if (currentHeader != null) {
                        String[] fields = line.split("\\s+");
                        if (fields.length == 7) {
                            int refStart = Integer.valueOf(fields[0]);
                            int refEnd = Integer.valueOf(fields[1]);
                            int altStart = Integer.valueOf(fields[2]);
                            int altEnd = Integer.valueOf(fields[3]);

                            Interval refInterval = refStart < refEnd ? new Interval(refContig, refStart, refEnd) : new Interval(refContig, refEnd, refStart);
                            Interval altInterval = altStart < altEnd ? new Interval(altContig, altStart, altEnd) : new Interval(altContig, altEnd, altStart);

                            if (reverse) {
                                imap.put(refInterval, altInterval);
                            } else {
                                imap.put(altInterval, refInterval);
                            }
                        }
                    }
                }
            }
        }

        return imap;
    }

    @Override
    public void execute() {
        log.info("Load alt vs. ref deltas...");
        IntervalTreeMap<Interval> imap = loadDeltas(false);

        Map<String, ContigInfo> cis = new HashMap<String, ContigInfo>();
        String sampleName = null, accession = null;

        for (SAMFileReader bamc : BAMSC) {
            if (sampleName == null || accession == null) {
                sampleName = bamc.getFileHeader().getReadGroups().iterator().next().getSample();
                accession = bamc.getFileHeader().getReadGroups().iterator().next().getReadGroupId();
            }

            for (SAMRecord read : bamc) {
                ContigInfo ci = cis.containsKey(accession + "." + read.getReadName()) ? cis.get(accession + "." + read.getReadName()) : new ContigInfo();

                ci.canonicalLocus = new Interval(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd());
                ci.numAlignmentsInCanonical++;
                ci.alignedToCanonical = read.getAlignmentStart() > 0;
                for (CigarElement ce : read.getCigar().getCigarElements()) {
                    ci.clippedInCanonical = ci.clippedInCanonical || ce.getOperator().equals(CigarOperator.H) || ce.getOperator().equals(CigarOperator.S);
                }

                if (read.getCigar().getCigarElements().size() == 1 && read.getCigar().getCigarElement(0).getOperator().equals(CigarOperator.M)) {
                    String md = read.getStringAttribute("MD");
                    try {
                        int length = Integer.parseInt(md);
                        if (length == read.getReadLength()) {
                            ci.perfectAlignmentInCanonical = true;
                        }
                    } catch (NumberFormatException e) {}
                }

                cis.put(accession + "." + read.getReadName(), ci);
            }
        }

        for (SAMFileReader bam0 : BAMS0) {
            for (SAMRecord read : bam0) {
                ContigInfo ci = cis.containsKey(accession + "." + read.getReadName()) ? cis.get(accession + "." + read.getReadName()) : new ContigInfo();

                ci.ref0Locus = new Interval(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd());
                ci.ref0ToCanonicalBlock = ci.ref0Locus;
                ci.ref0ToCanonicalExact = ci.ref0Locus;

                if (imap.getOverlapping(ci.ref0Locus).size() == 1) {
                    ci.ref0ToCanonicalBlock = imap.getOverlapping(ci.ref0Locus).iterator().next();

                    if (ci.ref0ToCanonicalBlock != null && imap.containsKey(ci.ref0ToCanonicalBlock)) {
                        Interval revLocus = imap.get(ci.ref0ToCanonicalBlock);
                        ci.ref0ToCanonicalExact = new Interval(ci.ref0ToCanonicalBlock.getSequence(),
                                                               ci.ref0ToCanonicalBlock.getStart() + ci.ref0ToCanonicalBlock.getStart() - revLocus.getStart(),
                                                               ci.ref0ToCanonicalBlock.getStart() + ci.ref0ToCanonicalBlock.getEnd() - revLocus.getStart());
                    }
                } else if (imap.getOverlapping(ci.ref0Locus).size() > 1) {
                    ci.ref0ToCanonicalBlock = null;
                    ci.ref0ToCanonicalExact = null;
                }

                ci.numAlignmentsInRef0++;
                ci.alignedToRef0 = read.getAlignmentStart() > 0;
                for (CigarElement ce : read.getCigar().getCigarElements()) {
                    ci.clippedInRef0 = ci.clippedInRef0 || ce.getOperator().equals(CigarOperator.H) || ce.getOperator().equals(CigarOperator.S);
                }

                if (read.getCigar().getCigarElements().size() == 1 && read.getCigar().getCigarElement(0).getOperator().equals(CigarOperator.M)) {
                    String md = read.getStringAttribute("MD");
                    try {
                        int length = Integer.parseInt(md);
                        if (length == read.getReadLength()) {
                            ci.perfectAlignmentInRef0 = true;
                        }
                    } catch (NumberFormatException e) {}
                }

                cis.put(accession + "." + read.getReadName(), ci);
            }
        }

        for (SAMFileReader bam1 : BAMS1) {
            for (SAMRecord read : bam1) {
                ContigInfo ci = cis.containsKey(accession + "." + read.getReadName()) ? cis.get(accession + "." + read.getReadName()) : new ContigInfo();

                ci.ref1Locus = new Interval(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd());
                ci.ref1ToCanonicalBlock = ci.ref1Locus;

                if (imap.getOverlapping(ci.ref1Locus).size() == 1) {
                    ci.ref1ToCanonicalBlock = imap.getOverlapping(ci.ref1Locus).iterator().next();
                } else if (imap.getOverlapping(ci.ref1Locus).size() > 1) {
                    ci.ref1ToCanonicalBlock = null;
                }

                if (ci.ref1ToCanonicalBlock != null && imap.containsKey(ci.ref1ToCanonicalBlock)) {
                    Interval revLocus = imap.get(ci.ref1ToCanonicalBlock);
                    ci.ref1ToCanonicalExact = new Interval(ci.ref1ToCanonicalBlock.getSequence(),
                                                           ci.ref1ToCanonicalBlock.getStart() + ci.ref1Locus.getStart() - revLocus.getStart(),
                                                           ci.ref1ToCanonicalBlock.getStart() + ci.ref1Locus.getEnd() - revLocus.getStart());
                }

                ci.sameChromosome = (ci.ref0ToCanonicalBlock != null && ci.ref1ToCanonicalBlock != null && ci.ref0ToCanonicalBlock.getSequence().equals(ci.ref1ToCanonicalBlock.getSequence()));
                ci.sameEffectiveLocus = (ci.ref0ToCanonicalBlock != null && ci.ref1ToCanonicalBlock != null && ci.ref0ToCanonicalBlock.intersects(ci.ref1ToCanonicalBlock));

                ci.numAlignmentsInRef1++;
                ci.alignedToRef1 = read.getAlignmentStart() > 0;
                for (CigarElement ce : read.getCigar().getCigarElements()) {
                    ci.clippedInRef1 = (ce.getOperator().equals(CigarOperator.H) || ce.getOperator().equals(CigarOperator.S));
                }

                if (read.getCigar().getCigarElements().size() == 1 && read.getCigar().getCigarElement(0).getOperator().equals(CigarOperator.M)) {
                    String md = read.getStringAttribute("MD");
                    try {
                        int length = Integer.parseInt(md);
                        if (length == read.getReadLength()) {
                            ci.perfectAlignmentInRef1 = true;
                        }
                    } catch (NumberFormatException e) {}
                }

                cis.put(accession + "." + read.getReadName(), ci);
            }
        }

        TableWriter tw = new TableWriter(out);

        log.info("Processing annotated contigs...");
        for (File ann : ANNS) {
            log.info("  {} {}", sampleName, accession);

            TableReader tr = new TableReader(ann);
            for (Map<String, String> te : tr) {
                int baseLength = te.get("seq").length();
                int kmerLength = te.get("kmerOrigin").length();

                int ref0 = 0;
                int ref1 = 0;
                int refBoth = 0;
                int refNone = 0;
                int refAmb = 0;

                for (int i = 0; i < kmerLength; i++) {
                    switch (te.get("kmerOrigin").charAt(i)) {
                        case '0': ref0++;    break;
                        case '1': ref1++;    break;
                        case 'B': refBoth++; break;
                        case '.': refNone++; break;
                        case 'A': refAmb++;  break;
                        default :            break;
                    }
                }

                boolean firstZeroSeen = false;
                boolean hasDiscontiguities = false;
                int numDiscontiguities = 0;

                for (int i = 1; i < kmerLength; i++) {
                    char kmerOrigin = te.get("kmerOrigin").charAt(i);
                    char contiguity = te.get("kmerContiguity").charAt(i);

                    if (kmerOrigin != '.') {
                        if (contiguity == '0') {
                            if (!firstZeroSeen) { firstZeroSeen = true; }
                            else {
                                hasDiscontiguities = true;
                                numDiscontiguities++;
                            }
                        }
                    }
                }

                int lrun0 = longestRun(te.get("kmerOrigin"), '0');
                int lrun1 = longestRun(te.get("kmerOrigin"), '1');
                int lrunBoth = longestRun(te.get("kmerOrigin"), 'B');
                int lrunNone = longestRun(te.get("kmerOrigin"), '.');
                int lrunAmb = longestRun(te.get("kmerOrigin"), 'A');
                int numSwitches = numSwitches(te.get("kmerOrigin"));

                Map<String, String> entry = new LinkedHashMap<String, String>();
                entry.put("sampleName", sampleName);
                entry.put("accession", accession);
                entry.put("contigName", te.get("contigName"));
                entry.put("baseLength", String.valueOf(baseLength));
                entry.put("kmerLength", String.valueOf(kmerLength));
                entry.put("ref0", String.valueOf(ref0));
                entry.put("ref1", String.valueOf(ref1));
                entry.put("refBoth", String.valueOf(refBoth));
                entry.put("refNone", String.valueOf(refNone));
                entry.put("refAmb", String.valueOf(refAmb));
                entry.put("lrun0", String.valueOf(lrun0));
                entry.put("lrun1", String.valueOf(lrun1));
                entry.put("lrunBoth", String.valueOf(lrunBoth));
                entry.put("lrunNone", String.valueOf(lrunNone));
                entry.put("lrunAmb", String.valueOf(lrunAmb));
                entry.put("numSwitches", String.valueOf(numSwitches));
                entry.put("hasDiscontiguities", hasDiscontiguities ? "1" : "0");
                entry.put("numDiscontiguities", String.valueOf(numDiscontiguities));

                if (cis.containsKey(accession + "." + te.get("contigName"))) {
                    ContigInfo ci = cis.get(accession + "." + te.get("contigName"));

                    entry.put("canonicalLocus", ci.canonicalLocus.toString());
                    entry.put("ref0Locus", ci.ref0Locus.toString());
                    entry.put("ref1Locus", ci.ref1Locus.toString());
                    entry.put("ref0ToCanonicalBlock", ci.ref0ToCanonicalBlock == null ? "NA" : ci.ref0ToCanonicalBlock.toString());
                    entry.put("ref1ToCanonicalBlock", ci.ref1ToCanonicalBlock == null ? "NA" : ci.ref1ToCanonicalBlock.toString());
                    entry.put("ref0ToCanonicalExact", ci.ref0ToCanonicalExact == null ? "NA" : ci.ref0ToCanonicalExact.toString());
                    entry.put("ref1ToCanonicalExact", ci.ref1ToCanonicalExact == null ? "NA" : ci.ref1ToCanonicalExact.toString());
                    entry.put("sameChromosome", ci.sameChromosome ? "1" : "0");
                    entry.put("sameEffectiveLocus", ci.sameEffectiveLocus ? "1" : "0");

                    Set<String> genes0 = new TreeSet<String>();
                    Set<String> genes1 = new TreeSet<String>();

                    if (ci.ref0ToCanonicalExact != null) {
                        for (GFF3Record gr : GFF3.getType("gene", GFF.getOverlapping(ci.ref0ToCanonicalExact))) {
                            genes0.add(gr.getAttribute("ID"));
                        }
                    }

                    if (ci.ref1ToCanonicalExact != null) {
                        for (GFF3Record gr : GFF3.getType("gene", GFF.getOverlapping(ci.ref1ToCanonicalExact))) {
                            genes1.add(gr.getAttribute("ID"));
                        }
                    }

                    entry.put("genesRef0", genes0.size() == 0 ? "NA" : Joiner.on(",").join(genes0));
                    entry.put("genesRef1", genes1.size() == 0 ? "NA" : Joiner.on(",").join(genes1));

                    entry.put("alignedToCanonical", ci.alignedToCanonical ? "1" : "0");
                    entry.put("clippedInCanonical", ci.clippedInCanonical ? "1" : "0");
                    entry.put("numAlignmentsInCanonical", String.valueOf(ci.numAlignmentsInCanonical));
                    entry.put("perfectAlignmentInCanonical", ci.perfectAlignmentInCanonical ? "1" : "0");
                    entry.put("alignedToRef0", ci.alignedToRef0 ? "1" : "0");
                    entry.put("clippedInRef0", ci.clippedInRef0 ? "1" : "0");
                    entry.put("numAlignmentsInRef0", String.valueOf(ci.numAlignmentsInRef0));
                    entry.put("perfectAlignmentInRef0", ci.perfectAlignmentInRef0 ? "1" : "0");
                    entry.put("alignedToRef1", ci.alignedToRef1 ? "1" : "0");
                    entry.put("clippedInRef1", ci.clippedInRef1 ? "1" : "0");
                    entry.put("numAlignmentsInRef1", String.valueOf(ci.numAlignmentsInRef1));
                    entry.put("perfectAlignmentInRef1", ci.perfectAlignmentInRef1 ? "1" : "0");
                } else {
                    entry.put("canonicalLocus", "NA");
                    entry.put("ref0Locus", "NA");
                    entry.put("ref1Locus", "NA");
                    entry.put("ref0ToCanonicalBlock", "NA");
                    entry.put("ref1ToCanonicalBlock", "NA");
                    entry.put("ref0ToCanonicalExact", "NA");
                    entry.put("ref1ToCanonicalExact", "NA");
                    entry.put("sameChromosome", "NA");
                    entry.put("sameEffectiveLocus", "NA");
                    entry.put("genesRef0", "NA");
                    entry.put("genesRef1", "NA");
                    entry.put("alignedToCanonical", "NA");
                    entry.put("clippedInCanonical", "NA");
                    entry.put("numAlignmentsInCanonical", "NA");
                    entry.put("perfectAlignmentInCanonical", "NA");
                    entry.put("alignedToRef0", "NA");
                    entry.put("clippedInRef0", "NA");
                    entry.put("numAlignmentsInRef0", "NA");
                    entry.put("perfectAlignmentInRef0", "NA");
                    entry.put("alignedToRef1", "NA");
                    entry.put("clippedInRef1", "NA");
                    entry.put("numAlignmentsInRef1", "NA");
                    entry.put("perfectAlignmentInRef1", "NA");
                }

                //entry.put("seq", te.get("seq"));
                //entry.put("kmerOrigin", te.get("kmerOrigin"));
                //entry.put("kmerContiguity", te.get("kmerContiguity"));
                //entry.put("kmerCoverage", te.get("kmerCoverage"));

                tw.addEntry(entry);
            }
        }
    }
}
