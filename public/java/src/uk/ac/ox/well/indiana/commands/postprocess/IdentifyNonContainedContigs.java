package uk.ac.ox.well.indiana.commands.postprocess;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;

import java.io.File;
import java.io.PrintStream;
import java.util.*;
import java.util.function.Predicate;

public class IdentifyNonContainedContigs extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAMs")
    public ArrayList<SamReader> BAMS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, SAMRecord> contigs = new HashMap<>();

        for (SamReader sam : BAMS) {
            for (SAMRecord sr : sam) {
                if (!contigs.containsKey(sr.getReadName())) {
                    contigs.put(sr.getReadName(), sr);
                } else {
                    SAMRecord old = contigs.get(sr.getReadName());

                    if (sr.getAlignmentEnd() - sr.getAlignmentStart() > old.getAlignmentEnd() - old.getAlignmentStart()) {
                        contigs.put(sr.getReadName(), sr);
                    }
                }
            }
        }

        IntervalTreeMap<String> t = new IntervalTreeMap<>();

        for (SAMRecord sr : contigs.values()) {
            Interval loc = new Interval(sr.getReferenceName(), sr.getAlignmentStart(), sr.getAlignmentEnd());
            t.put(loc, sr.getReadName());
        }

        List<String> containedReads = new ArrayList<>();

        for (SAMRecord sr : contigs.values()) {
            Interval loc = new Interval(sr.getReferenceName(), sr.getAlignmentStart(), sr.getAlignmentEnd());

            Collection<String> contained = t.getContained(loc);
            contained.remove(sr.getReadName());

            containedReads.addAll(contained);
        }

        for (String contigName : contigs.keySet()) {
            if (!containedReads.contains(contigName)) {
                out.println(contigName);
            }
        }
    }
}
