package uk.ac.ox.well.indiana.analyses.kmerSharing;

import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataFrame;
import uk.ac.ox.well.indiana.utils.io.utils.TableReader;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

public class ComputeRecoveryMatrix extends Tool {
    @Argument(fullName="referenceTable", shortName="rt", doc="Reference table (the table that contains *all* kmers)")
    public File REFERENCE_TABLE;

    @Argument(fullName="tables", shortName="t", doc="Tables of other samples' sequences")
    public ArrayList<File> TABLES;

    @Output
    public PrintStream out;

    @Override
    public int execute() {
        DataFrame<String, String, Integer> recovery = new DataFrame<String, String, Integer>(0);

        TableReader refTableReader = new TableReader(REFERENCE_TABLE);
        String refSampleName = REFERENCE_TABLE.getName().replaceAll("relatedSequences.", "").replaceAll(".table", "");
        for (HashMap<String, String> entry : refTableReader) {
            String geneName = entry.get("genes");

            recovery.set(refSampleName, geneName, recovery.get(refSampleName, geneName) + 1);
        }

        for (File table : TABLES) {
            if (!table.equals(REFERENCE_TABLE)) {
                TableReader tableReader = new TableReader(table);

                log.info("Loading '{}'", table.getName());

                String sampleName = table.getName().replaceAll("relatedSequences.", "").replaceAll(".table", "");
                for (HashMap<String, String> entry : tableReader) {
                    String geneName = entry.get("genes");

                    recovery.set(sampleName, geneName, recovery.get(sampleName, geneName) + 1);
                }
            }
        }

        DataFrame<String, String, Float> recoveryFraction = new DataFrame<String, String, Float>(0.0f);

        for (String sampleName : recovery.getRowNames()) {
            for (String geneName : recovery.getColNames()) {
                if (!sampleName.equalsIgnoreCase(refSampleName)) {
                    float fraction = ((float) recovery.get(sampleName, geneName)) / ((float) recovery.get(refSampleName, geneName));

                    recoveryFraction.set(sampleName, geneName, fraction);
                }
            }
        }

        out.println(recoveryFraction);

        return 0;
    }
}
