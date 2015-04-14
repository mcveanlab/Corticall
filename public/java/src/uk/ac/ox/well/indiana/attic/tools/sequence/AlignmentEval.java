package uk.ac.ox.well.indiana.attic.tools.sequence;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class AlignmentEval extends Module {
    @Argument(fullName="baseline", shortName="b", doc="Baseline BAM")
    public SAMFileReader BASELINE;

    @Argument(fullName="comparisons", shortName="c", doc="Comparison BAM(s)")
    public HashMap<String, SAMFileReader> COMPARISONS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, SAMRecordIterator> sis = new HashMap<String, SAMRecordIterator>();
        for (String key : COMPARISONS.keySet()) {
            SAMFileReader sfr = COMPARISONS.get(key);
            sis.put(key, sfr.iterator());
        }

        SAMRecordIterator bis = BASELINE.iterator();

        while (bis.hasNext()) {

        }
    }
}
