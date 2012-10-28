package uk.ac.ox.well.indiana.tools.cortex;

import com.google.common.base.Joiner;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.statistics.StatisticsOnStream;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class ViewCortexCoverage extends ViewCortexBase {
    private HashMap<String, HashMap<Integer, StatisticsOnStream>> stats = new HashMap<String, HashMap<Integer, StatisticsOnStream>>();

    @Override
    public int execute() {
        stats.put("total", new HashMap<Integer, StatisticsOnStream>());

        for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
            stats.get("total").put(color, new StatisticsOnStream());
        }

        log.info("Processing Cortex records");
        for (CortexRecord cr : CORTEX_GRAPH) {
            if (satisfiesConstraints(cr)) {
                String home = getKmerHomeContigName(cr);

                if (home == null || !"none".equals(home)) {
                    int[] coverages = cr.getCoverages();

                    for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
                        int coverage = coverages[color];

                        if (home != null) {
                            if (!stats.containsKey(home)) {
                                stats.put(home, new HashMap<Integer, StatisticsOnStream>());
                            }

                            if (!stats.get(home).containsKey(color)) {
                                stats.get(home).put(color, new StatisticsOnStream());
                            }

                            stats.get(home).get(color).push(coverage);
                        }

                        stats.get("total").get(color).push(coverage);
                    }
                }
            }
        }

        List<String> header = new ArrayList<String>();
        header.add("home");
        for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
            header.add(String.valueOf(color));
        }

        out.println(Joiner.on("\t").join(header));

        for (String home : stats.keySet()) {
            List<String> entry = new ArrayList<String>();

            entry.add(home);

            for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
                StatisticsOnStream st = stats.get(home).get(color);

                if (st == null) {
                    entry.add("0");
                } else {
                    entry.add(String.valueOf(st.getMean()));
                }
            }

            out.println(Joiner.on("\t").useForNull("0").join(entry));
        }

        return 0;
    }
}
