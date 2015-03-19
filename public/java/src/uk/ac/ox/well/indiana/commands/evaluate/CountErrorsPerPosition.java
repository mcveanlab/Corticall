package uk.ac.ox.well.indiana.commands.evaluate;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataTable;
import uk.ac.ox.well.indiana.utils.sequence.AlignmentUtils;

import java.io.PrintStream;
import java.util.Set;
import java.util.TreeSet;

public class CountErrorsPerPosition extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public SAMFileReader BAM;

    @Argument(fullName="maxNumReads", shortName="m", doc="Maximum number of reads to process (0 for all)")
    public Integer MAX_NUM_READS = 0;

    @Output
    public PrintStream out;

    private void initializeTable(DataTable t, int readLength) {
        for (boolean end : new boolean[] { false, true }) {
            for (int pos = 0; pos < readLength; pos++) {
                String pk = String.format("%b:%d", end, pos);
                t.set(pk, "isFirstEnd", end);
                t.set(pk, "position", pos);
                t.set(pk, "numMismatches", 0);
                t.set(pk, "numInsertions", 0);
                t.set(pk, "numDeletions", 0);
            }
        }
    }

    private Set<Integer> flip(Set<Integer> positions, int readLength) {
        Set<Integer> flippedPositions = new TreeSet<Integer>();

        for (int pos : positions) {
            flippedPositions.add(readLength - pos - 1);
        }

        return flippedPositions;
    }

    private void increment(DataTable t, boolean isFirstEnd, String column, Set<Integer> positions, SAMRecord read) {
        if (read.getReadNegativeStrandFlag()) {
            positions = flip(positions, read.getReadLength());
        }

        for (int pos : positions) {
            String pk = String.format("%b:%d", isFirstEnd, pos);
            int oldValue = t.has(pk, column) ? (Integer) t.get(pk, column) : 0;

            t.set(pk, "isFirstEnd", isFirstEnd);
            t.set(pk, "position", pos);
            t.set(pk, column, oldValue + 1);
        }
    }

    private void normalize(DataTable t, long numProcessedReads) {
        for (String pk : t.getPrimaryKeys()) {
            int numMismatches = t.has(pk, "numMismatches") ? (Integer) t.get(pk, "numMismatches") : 0;
            int numInsertions = t.has(pk, "numInsertions") ? (Integer) t.get(pk, "numInsertions") : 0;
            int numDeletions =  t.has(pk, "numDeletions")  ? (Integer) t.get(pk, "numDeletions")  : 0;

            t.set(pk, "rateMismatches", (double) numMismatches / (double) numProcessedReads);
            t.set(pk, "rateInsertions", (double) numInsertions / (double) numProcessedReads);
            t.set(pk, "rateDeletions",  (double) numDeletions  / (double) numProcessedReads);
        }
    }

    @Override
    public void execute() {
        DataTable t = new DataTable("errorPerPosition", "Table of errors per position in read");
        t.addColumns("isFirstEnd", "position", "numMismatches", "numInsertions", "numDeletions", "rateMismatches", "rateInsertions", "rateDeletions");

        boolean isInitialized = false;

        log.info("Processing reads...");
        long numReads = 0;
        long numProcessedReads = 0;
        for (SAMRecord read : BAM) {
            if (numReads % 1000000 == 0) {
                log.info("  {} reads", numReads);
            }
            numReads++;

            if (!read.getReadUnmappedFlag()) {
                numProcessedReads++;

                if (!isInitialized) {
                    initializeTable(t, read.getReadLength());
                    isInitialized = true;
                }

                boolean isFirstEnd = read.getFirstOfPairFlag();

                increment(t, isFirstEnd, "numMismatches", AlignmentUtils.getMismatchPositions(read), read);
                increment(t, isFirstEnd, "numInsertions", AlignmentUtils.getInsertionPositions(read), read);
                increment(t, isFirstEnd, "numDeletions", AlignmentUtils.getDeletionPositions(read), read);
            }

            if (MAX_NUM_READS > 0 && numReads > MAX_NUM_READS) {
                break;
            }
        }

        normalize(t, numProcessedReads);

        t.write(out);
    }
}
