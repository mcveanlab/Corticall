package uk.ac.ox.well.indiana.attic.analyses.nahr;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTreeMap;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;
import uk.ac.ox.well.indiana.utils.io.xmfa.XMFARecord;
import uk.ac.ox.well.indiana.utils.io.xmfa.XMFASequenceFile;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

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

    //@Argument(fullName="xmfa", shortName="xmfa", doc="XMFA file (Mauve alignment)")
    //public XMFASequenceFile XMFA;

    @Argument(fullName="deltas", shortName="d", doc="Delta")
    public File DELTA;

    @Argument(fullName="ref0", shortName="r0", doc="Ref 0")
    public String REF0;

    @Argument(fullName="ref1", shortName="r1", doc="Ref 1")
    public String REF1;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Output
    public PrintStream out;

    private class ContigInfo {
        public Interval ref0Locus;
        public Interval ref1Locus;
        public Interval ref1ToRef0Interval;
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

    /*
    private Map<CortexKmer, Set<String>> hashAlignments() {
        Map<CortexKmer, Set<String>> h = new HashMap<CortexKmer, Set<String>>();

        int processed = 0;
        for (XMFARecord xr : XMFA) {
            if (processed % (XMFA.getNumRecords() / 10) == 0) {
                log.info("  processed {}/{} (~{}%) records", processed, XMFA.getNumRecords(), String.format("%.2f", 100.0 * processed / XMFA.getNumRecords()));
            }
            processed++;

            String r0 = xr.containsKey(REF0) ? new String(xr.get(REF0).getBases()).replaceAll("-", "") : "";
            String r1 = xr.containsKey(REF1) ? new String(xr.get(REF1).getBases()) : "";

            for (int i = 0; i <= r0.length() - KMER_SIZE; i++) {
                CortexKmer kmer = new CortexKmer(r0.substring(i, i + KMER_SIZE));

                if (!h.containsKey(kmer)) {
                    h.put(kmer, new HashSet<String>());
                }

                h.get(kmer).add(r1);
            }
        }

        return h;
    }
    */

    private IntervalTreeMap<Interval> loadDeltas() {
        IntervalTreeMap<Interval> imap = new IntervalTreeMap<Interval>();

        LineReader lr = new LineReader(DELTA);

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

                        imap.put(altInterval, refInterval);
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
        IntervalTreeMap<Interval> imap = loadDeltas();

        Map<String, ContigInfo> cis = new HashMap<String, ContigInfo>();

        for (SAMFileReader bam1 : BAMS1) {
            String sampleName = bam1.getFileHeader().getReadGroups().iterator().next().getSample();

            for (SAMRecord read : bam1) {
                ContigInfo ci = cis.containsKey(sampleName + "." + read.getReadName()) ? cis.get(sampleName + "." + read.getReadName()) : new ContigInfo();

                ci.ref0Locus = new Interval(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd());
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

                ci.ref1Locus = new Interval(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd());

                log.info("{} {} {}", ci.ref1Locus, imap.get(ci.ref1Locus), imap.getOverlapping(ci.ref1Locus));
                for (Interval refLoci : imap.getOverlapping(ci.ref1Locus)) {
                    log.info("  {} {} {}", refLoci, ci.ref0Locus, ci.ref0Locus.intersects(refLoci));

                    if (ci.ref0Locus.intersects(refLoci)) {
                        ci.ref1ToRef0Interval = refLoci;
                        ci.sameEffectiveLocus = true;
                    }
                }

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
                if (te.get("contigName").equals("contig53046")) {
                    log.info("Hello!");
                }

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

                int numRef1KmersInSameStretchAsRef0Kmers = 0;
                int numRef1KmersNotInSameStretchAsRef0Kmers = 0;

                /*
                Set<String> relatedSeqs = new HashSet<String>();
                for (int i = 0; i < kmerLength; i++) {
                    char kmerOrigin = te.get("kmerOrigin").charAt(i);

                    if (kmerOrigin == '0' || kmerOrigin == 'B') {
                        CortexKmer kmer = new CortexKmer(te.get("seq").substring(i, i + KMER_SIZE));
                        if (h.containsKey(kmer)) { relatedSeqs.addAll(h.get(kmer)); }
                    }
                }

                for (int i = 0; i < kmerLength; i++) {
                    char kmerOrigin = te.get("kmerOrigin").charAt(i);

                    if (kmerOrigin == '1') {
                        String fw = te.get("seq").substring(i, i + KMER_SIZE);
                        String rc = SequenceUtils.reverseComplement(fw);

                        boolean found = false;

                        for (String relatedSeq : relatedSeqs) {
                            String seq = relatedSeq.replaceAll("-", "");
                            if (seq.contains(fw) || seq.contains(rc)) {
                                found = true;
                                break;
                            }
                        }

                        if (found) { numRef1KmersInSameStretchAsRef0Kmers++; }
                        else { numRef1KmersNotInSameStretchAsRef0Kmers++; }
                    }
                }
                */

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

                    entry.put("ref0Locus", ci.ref0Locus.toString());
                    entry.put("ref1Locus", ci.ref1Locus.toString());
                    entry.put("ref1ToRef0Interval", ci.ref1ToRef0Interval != null ? ci.ref1ToRef0Interval.toString() : "NA");
                    entry.put("sameEffectiveLocus", ci.sameEffectiveLocus ? "1" : "0");
                    entry.put("alignedToRef0", ci.alignedToRef0 ? "1" : "0");
                    entry.put("clippedInRef0", ci.clippedInRef0 ? "1" : "0");
                    entry.put("numAlignmentsInRef0", String.valueOf(ci.numAlignmentsInRef0));
                    entry.put("perfectAlignmentInRef0", ci.perfectAlignmentInRef0 ? "1" : "0");
                    entry.put("alignedToRef1", ci.alignedToRef1 ? "1" : "0");
                    entry.put("clippedInRef1", ci.clippedInRef1 ? "1" : "0");
                    entry.put("numAlignmentsInRef1", String.valueOf(ci.numAlignmentsInRef1));
                    entry.put("perfectAlignmentInRef1", ci.perfectAlignmentInRef1 ? "1" : "0");
                } else {
                    entry.put("ref0Locus", "NA");
                    entry.put("ref1Locus", "NA");
                    entry.put("ref1ToRef0Interval", "NA");
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

                entry.put("ref1KmersInRef0Stretches", String.valueOf(numRef1KmersInSameStretchAsRef0Kmers));
                entry.put("ref1KmersNotInRef0Stretches", String.valueOf(numRef1KmersNotInSameStretchAsRef0Kmers));

                entry.put("seq", te.get("seq"));
                entry.put("kmerOrigin", te.get("kmerOrigin"));
                entry.put("kmerContiguity", te.get("kmerContiguity"));

                tw.addEntry(entry);
            }
        }
    }
}
