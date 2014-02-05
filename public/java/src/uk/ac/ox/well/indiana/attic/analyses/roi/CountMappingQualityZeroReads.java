package uk.ac.ox.well.indiana.attic.analyses.roi;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.File;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

public class CountMappingQualityZeroReads extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public File BAM;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        SAMFileReader bamReader = new SAMFileReader(BAM);
        bamReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        Map<String, Integer> mq0 = new HashMap<String, Integer>();
        Map<String, Integer> unmapped = new HashMap<String, Integer>();
        Map<String, Integer> total = new HashMap<String, Integer>();

        for (SAMRecord record : bamReader) {
            String sample = record.getReadGroup().getSample();

            if (record.getMappingQuality() == 0) {
                if (!mq0.containsKey(sample)) {
                    mq0.put(sample, 0);
                }

                mq0.put(sample, mq0.get(sample) + 1);
            }

            if (record.getReadUnmappedFlag()) {
                if (!unmapped.containsKey(sample)) {
                    unmapped.put(sample, 0);
                }

                unmapped.put(sample, unmapped.get(sample) + 1);
            }

            if (!total.containsKey(sample)) {
                total.put(sample, 0);
            }

            total.put(sample, total.get(sample) + 1);
        }

        for (String sample : mq0.keySet()) {
            out.println(sample + "\t" + mq0.get(sample) + "\t" + unmapped.get(sample) + "\t" + total.get(sample));
        }
    }
}
