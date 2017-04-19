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

public class RemoveContainedContigs extends Module {
    @Argument(fullName="ref1", shortName="r1", doc="Reference 1")
    public IndexedFastaSequenceFile REF1;

    @Argument(fullName="bam1", shortName="b1", doc="BAM 1")
    public SamReader BAM1;

    @Argument(fullName="ref2", shortName="r2", doc="Reference 2")
    public IndexedFastaSequenceFile REF2;

    @Argument(fullName="bam2", shortName="b2", doc="BAM 2")
    public SamReader BAM2;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, SAMRecord> contigs = new HashMap<>();

        for (SamReader sam : Arrays.asList(BAM1, BAM2)) {
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

            containedReads.addAll(t.getContained(loc));
        }

        for (String contigName : contigs.keySet()) {
            if (!contigs.containsKey(contigName)) {
                out.println(contigName);
            }
        }
    }
}
