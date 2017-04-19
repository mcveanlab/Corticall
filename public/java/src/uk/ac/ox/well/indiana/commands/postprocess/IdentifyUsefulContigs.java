package uk.ac.ox.well.indiana.commands.postprocess;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.PrintStream;
import java.util.*;

public class IdentifyUsefulContigs extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAMs")
    public ArrayList<SamReader> BAMS;

    @Argument(fullName="minlength", shortName="l", doc="Minimum length to retain")
    public Integer MIN_LENGTH = 1000;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        List<String> containedReads = new ArrayList<>();
        Map<String, Integer> allContigs = new HashMap<>();

        for (SamReader sam : BAMS) {
            Map<String, SAMRecord> contigs = new HashMap<>();

            for (SAMRecord sr : sam) {
                allContigs.put(sr.getReadName(), sr.getReadLength());

                if (!contigs.containsKey(sr.getReadName())) {
                    contigs.put(sr.getReadName(), sr);
                } else {
                    SAMRecord old = contigs.get(sr.getReadName());

                    if (sr.getAlignmentEnd() - sr.getAlignmentStart() > old.getAlignmentEnd() - old.getAlignmentStart()) {
                        contigs.put(sr.getReadName(), sr);
                    }
                }
            }

            IntervalTreeMap<String> t = new IntervalTreeMap<>();

            for (SAMRecord sr : contigs.values()) {
                Interval loc = new Interval(sr.getReferenceName(), sr.getAlignmentStart(), sr.getAlignmentEnd());
                t.put(loc, sr.getReadName());
            }

            for (SAMRecord sr : contigs.values()) {
                Interval loc = new Interval(sr.getReferenceName(), sr.getAlignmentStart(), sr.getAlignmentEnd());

                Collection<String> contained = t.getContained(loc);
                contained.remove(sr.getReadName());

                containedReads.addAll(contained);
            }
        }

        for (String contigName : allContigs.keySet()) {
            if (!containedReads.contains(contigName) && allContigs.get(contigName) >= MIN_LENGTH) {
                out.println(contigName);
            }
        }
    }
}
