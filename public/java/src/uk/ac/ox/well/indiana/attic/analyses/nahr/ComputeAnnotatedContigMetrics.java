package uk.ac.ox.well.indiana.attic.analyses.nahr;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

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

    @Output
    public PrintStream out;

    private class ContigInfo {
        public boolean alignedToRef1 = false;
        public boolean alignedToRef2 = false;
        public boolean clippedInRef1 = false;
        public boolean clippedInRef2 = false;
        public int numAlignmentsInRef1 = 0;
        public int numAlignmentsInRef2 = 0;
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

    @Override
    public void execute() {
        Map<String, ContigInfo> cis = new HashMap<String, ContigInfo>();

        for (SAMFileReader bam1 : BAMS1) {
            for (SAMRecord read : bam1) {
                ContigInfo ci = cis.containsKey(read.getReadName()) ? cis.get(read.getReadName()) : new ContigInfo();

                ci.numAlignmentsInRef1++;
                ci.alignedToRef1 = read.getAlignmentStart() > 0;
                for (CigarElement ce : read.getCigar().getCigarElements()) {
                    ci.clippedInRef1 = (ce.getOperator().equals(CigarOperator.H) || ce.getOperator().equals(CigarOperator.S));
                }

                cis.put(read.getReadName(), ci);
            }
        }

        for (SAMFileReader bam2 : BAMS2) {
            for (SAMRecord read : bam2) {
                ContigInfo ci = cis.containsKey(read.getReadName()) ? cis.get(read.getReadName()) : new ContigInfo();

                ci.numAlignmentsInRef2++;
                ci.alignedToRef2 = read.getAlignmentStart() > 0;
                for (CigarElement ce : read.getCigar().getCigarElements()) {
                    ci.clippedInRef2 = (ce.getOperator().equals(CigarOperator.H) || ce.getOperator().equals(CigarOperator.S));
                }

                cis.put(read.getReadName(), ci);
            }
        }

        TableWriter tw = new TableWriter(out);

        log.info("Processing annotated contigs...");
        for (File ann : ANNS) {
            log.info("  {}", ann.getName());

            String sampleName = ann.getName().replaceAll(".contigs.unique.ann.fasta", "");

            TableReader tr = new TableReader(ann, new String[] {"contigName", "contig", "ann"});

            for (Map<String, String> te : tr) {
                int baseLength = te.get("contig").length();
                int kmerLength = te.get("ann").length();

                int ref1 = 0;
                int ref2 = 0;
                int refBoth = 0;
                int refNone = 0;

                for (int i = 0; i < kmerLength; i++) {
                    switch (te.get("ann").charAt(i)) {
                        case '0': ref1++;    break;
                        case '1': ref2++;    break;
                        case 'B': refBoth++; break;
                        case '.': refNone++; break;
                        default :            break;
                    }
                }

                int lrun1 = longestRun(te.get("ann"), '0');
                int lrun2 = longestRun(te.get("ann"), '1');
                int lrunBoth = longestRun(te.get("ann"), 'B');
                int lrunNone = longestRun(te.get("ann"), '.');
                int numSwitches = numSwitches(te.get("ann"));

                Map<String, String> entry = new LinkedHashMap<String, String>();
                entry.put("sampleName", sampleName);
                entry.put("contigName", te.get("contigName"));
                entry.put("baseLength", String.valueOf(baseLength));
                entry.put("kmerLength", String.valueOf(kmerLength));
                entry.put("ref1", String.valueOf(ref1));
                entry.put("ref2", String.valueOf(ref2));
                entry.put("refBoth", String.valueOf(refBoth));
                entry.put("refNone", String.valueOf(refNone));
                entry.put("lrun1", String.valueOf(lrun1));
                entry.put("lrun2", String.valueOf(lrun2));
                entry.put("lrunBoth", String.valueOf(lrunBoth));
                entry.put("lrunNone", String.valueOf(lrunNone));
                entry.put("numSwitches", String.valueOf(numSwitches));

                if (cis.containsKey(te.get("contigName"))) {
                    ContigInfo ci = cis.get(te.get("contigName"));

                    entry.put("alignedToRef1", ci.alignedToRef1 ? "1" : "0");
                    entry.put("clippedInRef1", ci.clippedInRef1 ? "1" : "0");
                    entry.put("numAlignmentsInRef1", String.valueOf(ci.numAlignmentsInRef1));
                    entry.put("alignedToRef2", ci.alignedToRef2 ? "1" : "0");
                    entry.put("clippedInRef2", ci.clippedInRef2 ? "1" : "0");
                    entry.put("numAlignmentsInRef2", String.valueOf(ci.numAlignmentsInRef2));
                } else {
                    entry.put("alignedToRef1", "0");
                    entry.put("clippedInRef1", "0");
                    entry.put("numAlignmentsInRef1", "0");
                    entry.put("alignedToRef2", "0");
                    entry.put("clippedInRef2", "0");
                    entry.put("numAlignmentsInRef2", "0");
                }

                tw.addEntry(entry);
            }
        }
    }
}
