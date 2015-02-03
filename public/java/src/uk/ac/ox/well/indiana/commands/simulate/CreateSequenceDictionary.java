package uk.ac.ox.well.indiana.commands.simulate;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;

import java.io.File;
import java.io.PrintStream;

public class CreateSequenceDictionary extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs")
    public File CONTIGS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        LineReader lr = new LineReader(CONTIGS);

        out.println("@HD\tVN:1.0\tSO:unsorted");

        String contigName = null;
        StringBuilder contigBuilder = new StringBuilder();
        while (lr.hasNext()) {
            String l = lr.getNextRecord();

            if (l.startsWith(">")) {
                if (contigName != null) {
                    out.println("@SQ\tSN:" + contigName + "\tLN:" + contigBuilder.length());
                }

                contigName = l.replaceAll(">", "").split("\\s+")[0];
                contigBuilder = new StringBuilder();
            } else {
                contigBuilder.append(l);
            }
        }

        out.println("@SQ\tSN:" + contigName + "\tLN:" + contigBuilder.length());

        //FastaSequenceIndex fsi = new FastaSequenceIndex(new File(CONTIGS.getAbsolutePath() + ".fai"));
        //log.info("fsi={}", fsi.size());
    }
}
