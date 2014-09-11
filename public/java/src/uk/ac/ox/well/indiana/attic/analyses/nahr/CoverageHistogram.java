package uk.ac.ox.well.indiana.attic.analyses.nahr;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.math.MoreMathUtils;

import java.io.PrintStream;

public class CoverageHistogram extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public SAMFileReader BAM;

    @Argument(fullName="binSize", shortName="s", doc="Bin size (bp)")
    public Integer BIN_SIZE = 2500;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        SAMSequenceDictionary dict = BAM.getFileHeader().getSequenceDictionary();

        IntervalTreeMap<Integer> hist = new IntervalTreeMap<Integer>();
        int unaligned = 0;

        for (SAMSequenceRecord sr : dict.getSequences()) {
            for (int binStart = 0; binStart < sr.getSequenceLength(); binStart += BIN_SIZE) {
                int binEnd = (binStart + BIN_SIZE) > sr.getSequenceLength() ? sr.getSequenceLength() : binStart + BIN_SIZE;

                Interval binInterval = new Interval(sr.getSequenceName(), binStart + 1, binEnd);

                hist.put(binInterval, 0);
            }
        }

        for (SAMRecord read : BAM) {
            Interval readInterval = new Interval(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentStart());

            for (int i = read.getAlignmentStart(); i < read.getAlignmentEnd(); i++) {
                Interval binInterval = new Interval(readInterval.getSequence(), MoreMathUtils.roundDown(i, BIN_SIZE) + 1, MoreMathUtils.roundUp(i, BIN_SIZE));

                if (hist.containsKey(binInterval)) {
                    hist.put(binInterval, hist.get(binInterval) + 1);
                } else {
                    unaligned++;
                }
            }

        }

        for (Interval binInterval : hist.keySet()) {
            out.printf("%s %d %d %.2f\n", binInterval.getSequence(), binInterval.getStart(), binInterval.getEnd(), (double) hist.get(binInterval) / (double) BIN_SIZE);
        }
        out.printf("NA 0 0 %.2f\n", (double) unaligned / (double) BIN_SIZE);
    }
}
