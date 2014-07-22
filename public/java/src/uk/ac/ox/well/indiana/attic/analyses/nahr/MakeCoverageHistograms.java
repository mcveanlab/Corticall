package uk.ac.ox.well.indiana.attic.analyses.nahr;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class MakeCoverageHistograms extends Module {
    @Argument(fullName="cov", shortName="c", doc="Coverage distribution")
    public HashMap<String, File> COV;

    @Argument(fullName="minCoverage", shortName="min", doc="Min coverage", required=false)
    public Integer MIN_COVERAGE = 0;

    @Argument(fullName="maxCoverage", shortName="max", doc="Max coverage", required=false)
    public Integer MAX_COVERAGE = 400;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, Map<Integer, Integer>> covHist = new HashMap<String, Map<Integer, Integer>>();

        //int minCoverage = 100;
        //int maxCoverage = 0;

        for (String sample : COV.keySet()) {
            File covDist = COV.get(sample);

            covHist.put(sample, new TreeMap<Integer, Integer>());

            TableReader tr = new TableReader(covDist);

            for (Map<String, String> te : tr) {
                int coverage = Integer.valueOf(te.get("sample"));

                //if (coverage < minCoverage) { minCoverage = coverage; }
                //if (coverage > maxCoverage) { maxCoverage = coverage; }

                if (!covHist.get(sample).containsKey(coverage)) {
                    covHist.get(sample).put(coverage, 1);
                } else {
                    covHist.get(sample).put(coverage, covHist.get(sample).get(coverage) + 1);
                }
            }
        }

        TableWriter tw = new TableWriter(out);

        for (int coverage = MIN_COVERAGE; coverage <= MAX_COVERAGE; coverage++) {
            Map<String, String> row = new LinkedHashMap<String, String>();

            row.put("coverage", String.valueOf(coverage));

            for (String sample : COV.keySet()) {
                int count = covHist.get(sample).containsKey(coverage) ? covHist.get(sample).get(coverage) : 0;

                row.put(sample, String.valueOf(count));
            }

            tw.addEntry(row);
        }
    }
}
