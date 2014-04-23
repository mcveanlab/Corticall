package uk.ac.ox.well.indiana.attic.analyses.nahr;

import com.google.common.base.Joiner;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTreeMap;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
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

    @Argument(fullName="bams1", shortName="b1", doc="BAMs")
    public ArrayList<SAMFileReader> BAMS1;

    @Argument(fullName="bams2", shortName="b2", doc="BAMs")
    public ArrayList<SAMFileReader> BAMS2;

    @Argument(fullName="deltas", shortName="d", doc="Delta")
    public ArrayList<File> DELTAS;

    @Argument(fullName="gff", shortName="g", doc="GFF file")
    public GFF3 GFF;

    @Output
    public PrintStream out;

    private class ContigInfo {
        public Interval originalRef0Locus;
        public Interval originalRef1Locus;
        public Interval convertedRef0Locus;
        public Interval convertedRef1Locus;
        public Interval exactRef0Locus;
        public Interval exactRef1Locus;
        public boolean sameChromosome = false;
        public boolean sameEffectiveLocus = false;
        public boolean alignedToRef0 = false;
        public boolean alignedToRef1 = false;
        public boolean clippedInRef0 = false;
        public boolean clippedInRef1 = false;
        public boolean perfectAlignmentInRef0 = false;
        public boolean perfectAlignmentInRef1 = false;
        public int numAlignmentsInRef0 = 0;
        public int numAlignmentsInRef1 = 0;
        public boolean mapsToSimilarLoci = false;
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
        //log.info("Loading XMFA alignment file...");
        //Map<CortexKmer, Set<String>> h = hashAlignments();

        log.info("Load alt vs. ref deltas...");
        IntervalTreeMap<Interval> imap = loadDeltas(false);
        IntervalTreeMap<Interval> rmap = loadDeltas(true);

        Map<String, ContigInfo> cis = new HashMap<String, ContigInfo>();

        for (SAMFileReader bam1 : BAMS1) {
            String sampleName = bam1.getFileHeader().getReadGroups().iterator().next().getSample();

            for (SAMRecord read : bam1) {
                ContigInfo ci = cis.containsKey(sampleName + "." + read.getReadName()) ? cis.get(sampleName + "." + read.getReadName()) : new ContigInfo();

                ci.originalRef0Locus = new Interval(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd());
                ci.convertedRef0Locus = ci.originalRef0Locus;
                ci.exactRef0Locus = ci.originalRef0Locus;

                if (imap.getOverlapping(ci.originalRef0Locus).size() == 1) {
                    ci.convertedRef0Locus = imap.getOverlapping(ci.originalRef0Locus).iterator().next();

                    if (ci.convertedRef0Locus != null && rmap.containsKey(ci.convertedRef0Locus)) {
                        Interval revLocus = rmap.get(ci.convertedRef0Locus);
                        ci.exactRef0Locus = new Interval(ci.convertedRef0Locus.getSequence(),
                                                         ci.convertedRef0Locus.getStart() + ci.originalRef0Locus.getStart() - revLocus.getStart(),
                                                         ci.convertedRef0Locus.getStart() + ci.originalRef0Locus.getEnd() - revLocus.getStart());
                    }
                } else if (imap.getOverlapping(ci.originalRef0Locus).size() > 1) {
                    ci.convertedRef0Locus = null;
                    ci.exactRef0Locus = null;
                }

                ci.numAlignmentsInRef0++;
                ci.alignedToRef0 = read.getAlignmentStart() > 0;
                for (CigarElement ce : read.getCigar().getCigarElements()) {
                    ci.clippedInRef0 = (ce.getOperator().equals(CigarOperator.H) || ce.getOperator().equals(CigarOperator.S));
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

                cis.put(sampleName + "." + read.getReadName(), ci);
            }
        }

        for (SAMFileReader bam2 : BAMS2) {
            String sampleName = bam2.getFileHeader().getReadGroups().iterator().next().getSample();

            for (SAMRecord read : bam2) {
                ContigInfo ci = cis.containsKey(sampleName + "." + read.getReadName()) ? cis.get(sampleName + "." + read.getReadName()) : new ContigInfo();

                ci.originalRef1Locus = new Interval(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd());
                ci.convertedRef1Locus = ci.originalRef1Locus;

                if (imap.getOverlapping(ci.originalRef1Locus).size() == 1) {
                    ci.convertedRef1Locus = imap.getOverlapping(ci.originalRef1Locus).iterator().next();
                } else if (imap.getOverlapping(ci.originalRef1Locus).size() > 1) {
                    ci.convertedRef1Locus = null;
                }

                //I [2014-04-23 02:07 19499] Supercontig_1.1:1742-1983 Pf3D7_14_v3:2650339-2705858 Supercontig_1.1:1524-57055

                if (ci.convertedRef1Locus != null && rmap.containsKey(ci.convertedRef1Locus)) {
                    Interval revLocus = rmap.get(ci.convertedRef1Locus);
                    ci.exactRef1Locus = new Interval(ci.convertedRef1Locus.getSequence(),
                                                     ci.convertedRef1Locus.getStart() + ci.originalRef1Locus.getStart() - revLocus.getStart(),
                                                     ci.convertedRef1Locus.getStart() + ci.originalRef1Locus.getEnd() - revLocus.getStart());
                }

                ci.sameChromosome = (ci.convertedRef0Locus != null && ci.convertedRef1Locus != null && ci.convertedRef0Locus.getSequence().equals(ci.convertedRef1Locus.getSequence()));
                ci.sameEffectiveLocus = (ci.convertedRef0Locus != null && ci.convertedRef1Locus != null && ci.convertedRef0Locus.intersects(ci.convertedRef1Locus));

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

                cis.put(sampleName + "." + read.getReadName(), ci);
            }
        }

        TableWriter tw = new TableWriter(out);

        log.info("Processing annotated contigs...");
        for (File ann : ANNS) {
            String sampleName = ann.getName().replaceAll(".contigs.unique.ann", "");

            log.info("  {}", sampleName);

            TableReader tr = new TableReader(ann);
            for (Map<String, String> te : tr) {
                int baseLength = te.get("seq").length();
                int kmerLength = te.get("kmerOrigin").length();

                int ref0 = 0;
                int ref1 = 0;
                int refBoth = 0;
                int refNone = 0;

                for (int i = 0; i < kmerLength; i++) {
                    switch (te.get("kmerOrigin").charAt(i)) {
                        case '0': ref0++;    break;
                        case '1': ref1++;    break;
                        case 'B': refBoth++; break;
                        case '.': refNone++; break;
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
                int numSwitches = numSwitches(te.get("kmerOrigin"));

                Map<String, String> entry = new LinkedHashMap<String, String>();
                entry.put("sampleName", sampleName);
                entry.put("contigName", te.get("contigName"));
                entry.put("baseLength", String.valueOf(baseLength));
                entry.put("kmerLength", String.valueOf(kmerLength));
                entry.put("ref0", String.valueOf(ref0));
                entry.put("ref1", String.valueOf(ref1));
                entry.put("refBoth", String.valueOf(refBoth));
                entry.put("refNone", String.valueOf(refNone));
                entry.put("lrun0", String.valueOf(lrun0));
                entry.put("lrun1", String.valueOf(lrun1));
                entry.put("lrunBoth", String.valueOf(lrunBoth));
                entry.put("lrunNone", String.valueOf(lrunNone));
                entry.put("numSwitches", String.valueOf(numSwitches));
                entry.put("hasDiscontiguities", hasDiscontiguities ? "1" : "0");
                entry.put("numDiscontiguities", String.valueOf(numDiscontiguities));

                if (cis.containsKey(sampleName + "." + te.get("contigName"))) {
                    ContigInfo ci = cis.get(sampleName + "." + te.get("contigName"));

                    entry.put("originalRef0Locus", ci.originalRef0Locus.toString());
                    entry.put("originalRef1Locus", ci.originalRef1Locus.toString());
                    entry.put("convertedRef0Locus", ci.convertedRef0Locus == null ? "NA" : ci.convertedRef0Locus.toString());
                    entry.put("convertedRef1Locus", ci.convertedRef1Locus == null ? "NA" : ci.convertedRef1Locus.toString());
                    entry.put("exactRef0Locus", ci.exactRef0Locus == null ? "NA" : ci.exactRef0Locus.toString());
                    entry.put("exactRef1Locus", ci.exactRef1Locus == null ? "NA" : ci.exactRef1Locus.toString());
                    entry.put("sameChromosome", ci.sameChromosome ? "1" : "0");
                    entry.put("sameEffectiveLocus", ci.sameEffectiveLocus ? "1" : "0");

                    Set<String> genes0 = new TreeSet<String>();
                    Set<String> genes1 = new TreeSet<String>();

                    if (ci.exactRef0Locus != null) {
                        for (GFF3Record gr : GFF3.getType("gene", GFF.getOverlapping(ci.exactRef0Locus))) {
                            genes0.add(gr.getAttribute("ID"));
                        }
                    }

                    if (ci.exactRef1Locus != null) {
                        for (GFF3Record gr : GFF3.getType("gene", GFF.getOverlapping(ci.exactRef1Locus))) {
                            genes1.add(gr.getAttribute("ID"));
                        }
                    }

                    entry.put("genesRef0", genes0.size() > 0 ? "NA" : Joiner.on(",").join(genes0));
                    entry.put("genesRef1", genes1.size() > 0 ? "NA" : Joiner.on(",").join(genes1));

//                    log.info("g1: {}", Joiner.on(",").join(genes0));
//                    log.info("g2: {}", Joiner.on(",").join(genes1));

                    entry.put("alignedToRef0", ci.alignedToRef0 ? "1" : "0");
                    entry.put("clippedInRef0", ci.clippedInRef0 ? "1" : "0");
                    entry.put("numAlignmentsInRef0", String.valueOf(ci.numAlignmentsInRef0));
                    entry.put("perfectAlignmentInRef0", ci.perfectAlignmentInRef0 ? "1" : "0");
                    entry.put("alignedToRef1", ci.alignedToRef1 ? "1" : "0");
                    entry.put("clippedInRef1", ci.clippedInRef1 ? "1" : "0");
                    entry.put("numAlignmentsInRef1", String.valueOf(ci.numAlignmentsInRef1));
                    entry.put("perfectAlignmentInRef1", ci.perfectAlignmentInRef1 ? "1" : "0");
                } else {
                    entry.put("originalRef0Locus", "NA");
                    entry.put("originalRef1Locus", "NA");
                    entry.put("convertedRef0Locus", "NA");
                    entry.put("convertedRef1Locus", "NA");
                    entry.put("exactRef0Locus", "NA");
                    entry.put("exactRef1Locus", "NA");
                    entry.put("sameChromosome", "NA");
                    entry.put("sameEffectiveLocus", "NA");
                    entry.put("alignedToRef0", "NA");
                    entry.put("clippedInRef0", "NA");
                    entry.put("numAlignmentsInRef0", "NA");
                    entry.put("perfectAlignmentInRef0", "NA");
                    entry.put("alignedToRef1", "NA");
                    entry.put("clippedInRef1", "NA");
                    entry.put("numAlignmentsInRef1", "NA");
                    entry.put("perfectAlignmentInRef1", "NA");
                }

                entry.put("seq", te.get("seq"));
                entry.put("kmerOrigin", te.get("kmerOrigin"));
                entry.put("kmerContiguity", te.get("kmerContiguity"));

                tw.addEntry(entry);
            }
        }
    }
}
