package uk.ac.ox.well.indiana.commands.simulate;

import be.ac.ulg.montefiore.run.jahmm.*;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchLearner;
import be.ac.ulg.montefiore.run.jahmm.toolbox.KullbackLeiblerDistanceCalculator;
import be.ac.ulg.montefiore.run.jahmm.toolbox.MarkovGenerator;
import com.google.common.base.Joiner;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.sequence.AlignmentUtils;

import java.io.PrintStream;
import java.util.*;

public class ComputeBaseAndFragmentErrorRates extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public SAMFileReader BAM;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TableWriter tw = new TableWriter(out);

        log.info("Processing reads...");
        for (SAMSequenceRecord ssr : BAM.getFileHeader().getSequenceDictionary().getSequences()) {
            log.info("  {}", ssr.getSequenceName());

            int[] numReadsWithErrors = new int[ssr.getSequenceLength()];
            int[] numReads = new int[ssr.getSequenceLength()];
            int[] numFragmentsWithErrors = new int[ssr.getSequenceLength()];
            int[] numFragments = new int[ssr.getSequenceLength()];

            SAMRecordIterator sri = BAM.queryOverlapping(ssr.getSequenceName(), 0, ssr.getSequenceLength());
            Map<String, SAMRecord> seenReads = new HashMap<String, SAMRecord>();

            while (sri.hasNext()) {
                SAMRecord sr = sri.next();

                Set<Integer> errorPositions = AlignmentUtils.getDifferencesPositions(sr);

                for (int pos = 0; pos < sr.getReadLength(); pos++) {
                    int refpos = sr.getReferencePositionAtReadPosition(pos);

                    if (refpos > 0 && refpos < numReads.length) {
                        if (errorPositions.contains(pos)) {
                            numReadsWithErrors[refpos - 1]++;
                        }
                        numReads[refpos - 1]++;
                    }
                }

                if (!seenReads.containsKey(sr.getReadName())) {
                    seenReads.put(sr.getReadName(), sr);
                } else {
                    SAMRecord srm = seenReads.get(sr.getReadName());

                    Set<Integer> errorPositionsM = AlignmentUtils.getDifferencesPositions(srm);

                    if (errorPositions.size() + errorPositionsM.size() > 0) {
                        numFragmentsWithErrors[srm.getAlignmentStart() - 1]++;
                    }
                    numFragments[srm.getAlignmentStart() - 1]++;

                    seenReads.remove(sr.getReadName());
                }
            }

            sri.close();

            for (int pos = 0; pos < numReads.length; pos++) {
                Map<String, String> te = new LinkedHashMap<String, String>();
                te.put("chr", ssr.getSequenceName());
                te.put("pos", String.valueOf(pos));
                te.put("numReadsWithErrors", String.valueOf(numReadsWithErrors[pos]));
                te.put("numReads", String.valueOf(numReads[pos]));
                te.put("numFragmentsWithErrors", String.valueOf(numFragmentsWithErrors[pos]));
                te.put("numFragments", String.valueOf(numFragments[pos]));

                tw.addEntry(te);
            }
        }
    }
}
