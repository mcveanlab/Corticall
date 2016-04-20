package uk.ac.ox.well.indiana.commands.evaluate;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataTable;

import java.io.PrintStream;
import java.util.*;

public class RealisticLanderWatermanStats extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public SAMFileReader BAM;

    @Argument(fullName="numReads", shortName="r", doc="Num reads")
    public ArrayList<Integer> NUM_READS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        log.info("Loading reads...");

        List<SAMRecord> reads = new ArrayList<SAMRecord>();

        for (SAMRecord sr : BAM) {
            reads.add(sr);
        }

        log.info("  {} reads loaded", reads.size());

        Collections.shuffle(reads);

        log.info("  shuffled");

        DataTable dt = new DataTable("LanderWaterman", "Realistic Lander-Waterman stats", "reads", "bp_used", "bp_total", "pct_genome");

        for (int numReads : NUM_READS) {
            log.info("Selecting {} reads...", numReads);

            Map<String, boolean[]> utilizationMask = new HashMap<String, boolean[]>();

            for (SAMSequenceRecord ssr : BAM.getFileHeader().getSequenceDictionary().getSequences()) {
                utilizationMask.put(ssr.getSequenceName(), new boolean[ssr.getSequenceLength()]);
            }

            for (int i = 0; i < numReads; i++) {
                SAMRecord read = reads.get(i);

                for (int j = read.getAlignmentStart(); j < read.getAlignmentEnd(); j++) {
                    utilizationMask.get(read.getReferenceName())[j] = true;
                }
            }

            int basesUsed = 0, basesTotal = 0;

            for (String refName : utilizationMask.keySet()) {
                for (int i = 0; i < utilizationMask.get(refName).length; i++) {
                    if (utilizationMask.get(refName)[i]) {
                        basesUsed++;
                    }

                    basesTotal++;
                }
            }

            dt.set("lw" + numReads, "reads", numReads);
            dt.set("lw" + numReads, "bp_used", basesUsed);
            dt.set("lw" + numReads, "bp_total", basesTotal);
            dt.set("lw" + numReads, "pct_genome", (float) basesUsed / (float) basesTotal);

            log.info("  reads: {}, bp_used: {}, bp_total: {}, pct_genome: {}", numReads, basesUsed, basesTotal, (float) basesUsed / (float) basesTotal);

        }

        out.println(dt);
    }
}
