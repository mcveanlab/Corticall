package uk.ac.ox.well.indiana.analyses.reconstruction;

import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataFrame;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

public class CombineSampleRelatednessMatrices extends Tool {
    @Argument(fullName="sampleRelatednessMatrices", shortName="srm", doc="Sample relatedness matrices")
    public ArrayList<File> MATRICES;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        DataFrame<String, String, Float> relatedness = new DataFrame<String, String, Float>(0.0f);

        for (File matrix : MATRICES) {
            DataFrame<String, String, Float> d = new DataFrame<String, String, Float>(matrix, 0.0f);

            relatedness.addMatrix(d);
        }

        out.println(relatedness);
    }
}
