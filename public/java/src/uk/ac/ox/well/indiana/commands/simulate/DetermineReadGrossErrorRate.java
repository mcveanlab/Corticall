package uk.ac.ox.well.indiana.commands.simulate;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

public class DetermineReadGrossErrorRate extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public ArrayList<SAMFileReader> BAMS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TableWriter tw = new TableWriter(out);

        for (SAMFileReader bam : BAMS) {
            long basesSeen = 0;
            long nmSeen = 0;

            for (SAMRecord read : bam) {
                if (read.getAttribute("NM") != null) {
                    int nm = read.getIntegerAttribute("NM");

                    basesSeen += read.getReadLength();
                    nmSeen += nm;
                }
            }

            Map<String, String> te = new LinkedHashMap<String, String>();
            te.put("sample", bam.getFileHeader().getReadGroups().get(0).getSample());
            te.put("accession", bam.getFileHeader().getReadGroups().get(0).getReadGroupId());
            te.put("basesSeen", String.valueOf(basesSeen));
            te.put("nmSeen", String.valueOf(nmSeen));
            te.put("errorRate", String.valueOf((float) nmSeen / (float) basesSeen));

            tw.addEntry(te);
        }
    }
}
