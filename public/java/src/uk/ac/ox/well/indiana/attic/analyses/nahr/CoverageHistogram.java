package uk.ac.ox.well.indiana.attic.analyses.nahr;

import htsjdk.samtools.*;
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

        IntervalTreeMap<Float> hist = new IntervalTreeMap<Float>();
        int unaligned = 0;

        log.info("Processing chromosomes...");
        for (SAMSequenceRecord ssr : dict.getSequences()) {
            log.info("  {} ({} bp)", ssr.getSequenceName(), ssr.getSequenceLength());

            for (int binStart = 0; binStart < ssr.getSequenceLength(); binStart += BIN_SIZE) {
                int binEnd = (binStart + BIN_SIZE) > ssr.getSequenceLength() ? ssr.getSequenceLength() : binStart + BIN_SIZE;

                Interval binInterval = new Interval(ssr.getSequenceName(), binStart + 1, binEnd);

                hist.put(binInterval, 0.0f);

                SAMRecordIterator sri = BAM.queryOverlapping(ssr.getSequenceName(), binStart + 1, binEnd);
                int covSum = 0;
                //int[] covs = new int[binEnd - (binStart + 1)];

                while (sri.hasNext()) {
                    SAMRecord sr = sri.next();

                    for (int i = sr.getAlignmentStart(); i < sr.getAlignmentEnd(); i++) {
                        //int binPos = binStart + 1 - i;

                        //if (binPos >= binStart + 1 && binPos < binEnd) {
                        if (i >= binStart + 1 && i < binEnd) {
                            //covs[binPos]++;

                            covSum++;
                        }
                    }
                }

                hist.put(binInterval, (float) covSum / (float) BIN_SIZE);

                sri.close();
            }
        }

        /*
        for (SAMRecord read : BAM) {
            Interval readInterval = new Interval(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentStart());

            int i = read.getAlignmentStart();
            Interval binInterval = new Interval(readInterval.getSequence(), MoreMathUtils.roundDown(i, BIN_SIZE) + 1, MoreMathUtils.roundUp(i, BIN_SIZE));

            if (hist.containsKey(binInterval)) {
                hist.put(binInterval, hist.get(binInterval) + 1);
            } else {
                unaligned++;
            }
        }
        */

        for (Interval binInterval : hist.keySet()) {
            out.printf("%s %d %d %.2f\n", binInterval.getSequence(), binInterval.getStart(), binInterval.getEnd(), hist.get(binInterval));
        }
        out.printf("NA 0 0 %d\n", unaligned);
    }
}
