package uk.ac.ox.well.indiana.attic.analyses.kmerSharing;

import com.google.common.base.Joiner;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataFrame;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.performance.PerformanceUtils;
import uk.ac.ox.well.indiana.utils.statistics.PCA;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Map;

public class CreateRelatednessMatrix extends Module {
    @Argument(fullName="tables", shortName="t", doc="Tables of other samples' sequences")
    public ArrayList<File> TABLES;

    @Argument(fullName="onGenes", shortName="og", doc="Create matrix on genes")
    public Boolean ON_GENES = false;

    @Output
    public PrintStream out;

    @Output(fullName="pcaOut", shortName="po", doc="The PCA output file")
    public PrintStream pcaOut;

    @Override
    public void execute() {
        DataFrame<String, String, Float> relatedness = new DataFrame<String, String, Float>(0.0f);

        for (File table : TABLES) {
            log.info("Reading table '{}'", table.getAbsolutePath());

            TableReader refTableReader = new TableReader(table);

            String colName = table.getName().replaceAll("relatedSequences.", "").replaceAll(".table", "");
            for (Map<String, String> entry : refTableReader) {
                if (ON_GENES) {
                    String geneName = entry.get("genes");

                    relatedness.set(entry.get("kmer"), colName + "." + geneName, 1.0f);
                } else {
                    relatedness.set(entry.get("kmer"), colName, 1.0f);
                }
            }
        }

        log.info("Writing relatedness matrix to disk");

        out.println("\t" + Joiner.on("\t").join(relatedness.getColNames()));

        for (String rowName : relatedness.getRowNames()) {
            ArrayList<String> entry = new ArrayList<String>();
            entry.add(rowName);

            for (String colName : relatedness.getColNames()) {
                entry.add(relatedness.get(rowName, colName).toString());
            }

            out.println(Joiner.on("\t").join(entry));
        }

        log.info("Performance: {}", PerformanceUtils.getMemoryUsageStats());

        log.info("Computing PCA");
        PCA<String, String> pca = new PCA<String, String>(relatedness, false, false);

        pcaOut.println(pca.getRotationMatrix());
        double[] sds = pca.getStandardDeviations();

        for (int i = 0; i < sds.length; i++) {
            pcaOut.println(i + " " + sds[i]);
        }
    }
}
