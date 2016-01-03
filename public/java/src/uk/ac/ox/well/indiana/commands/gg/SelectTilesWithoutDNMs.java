package uk.ac.ox.well.indiana.commands.gg;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class SelectTilesWithoutDNMs extends Module {
    @Argument(fullName="featureTable", shortName="f", doc="Feature table")
    public ArrayList<File> FEATURE_TABLES;

    @Argument(fullName="variantTable", shortName="v", doc="Variant table")
    public File VARIANT_TABLE;

    @Argument(fullName="limit", shortName="l", doc="Limit")
    public Integer LIMIT = 1000;

    @Argument(fullName="seed", shortName="s", doc="Seed")
    public Long SEED = System.currentTimeMillis();

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        log.info("Loading features...");
        IntervalTreeMap<Map<String, String>> features = new IntervalTreeMap<Map<String, String>>();

        for (File featureTable : FEATURE_TABLES) {
            TableReader tr = new TableReader(featureTable);

            for (Map<String, String> te : tr) {
                String chr = te.get("chr");
                int start = Integer.valueOf(te.get("start"));
                int stop = Integer.valueOf(te.get("stop"));

                Interval interval = new Interval(chr, start, stop);

                features.put(interval, te);
            }
        }

        log.info("  features: {}", features.values().size());

        log.info("Loading variant loci...");
        TableReader tr = new TableReader(VARIANT_TABLE);

        for (Map<String, String> ve : tr) {
            String locus = ve.get("locus");

            if (!locus.contains("NA")) {
                String[] pieces = locus.split("[:-]");
                String vchr = pieces[0];
                int vstart = Integer.valueOf(pieces[1]);
                int vstop = Integer.valueOf(pieces[2]);

                Interval vinterval = new Interval(vchr, vstart, vstop);

                if (features.containsOverlapping(vinterval)) {
                    for (Map<String, String> te : features.getOverlapping(vinterval)) {
                        String chr = te.get("chr");
                        int start = Integer.valueOf(te.get("start"));
                        int stop = Integer.valueOf(te.get("stop"));

                        Interval interval = new Interval(chr, start, stop);
                        features.remove(interval);
                    }
                }
            }
        }

        log.info("  features: {}", features.values().size());

        Random rng = new Random();
        rng.setSeed(SEED);

        Set<Integer> ids = new TreeSet<Integer>();
        int limit = features.values().size() < LIMIT ? features.values().size() : LIMIT;

        for (int i = 0; i < limit; i++) {
            int index;
            do {
                index = rng.nextInt(features.values().size());
            } while (ids.contains(index));

            ids.add(index);
        }

        TableWriter tw = new TableWriter(out);

        List<Map<String, String>> tes = new ArrayList<Map<String, String>>(features.values());
        for (Integer index : ids) {
            Map<String, String> te = tes.get(index);
            te.put("index", String.valueOf(index));

            tw.addEntry(te);
        }
    }
}
