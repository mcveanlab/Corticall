package uk.ac.ox.well.indiana.commands.gg;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class AnnotateCalls extends Module {
    @Argument(fullName="gffs", shortName="g", doc="GFFs")
    public ArrayList<GFF3> GFFS;

    @Argument(fullName="properties", shortName="p", doc="Sequence properties")
    public ArrayList<File> PROPS;

    @Argument(fullName="calls", shortName="c", doc="Calls")
    public File CALLS;

    @Output
    public PrintStream out;

    private IntervalTreeMap<Map<String, String>> loadProperties() {
        IntervalTreeMap<Map<String, String>> tm = new IntervalTreeMap<Map<String, String>>();

        for (File file : PROPS) {
            TableReader tr = new TableReader(file);

            for (Map<String, String> te : tr) {
                String chr = te.get("chr");
                int start = Integer.valueOf(te.get("start"));
                int stop = Integer.valueOf(te.get("stop"));

                Interval interval = new Interval(chr, start, stop);

                te.remove("chr");
                te.remove("start");
                te.remove("stop");

                tm.put(interval, te);
            }
        }

        return tm;
    }

    @Override
    public void execute() {
        IntervalTreeMap<Map<String, String>> p = loadProperties();

        Set<String> props = new TreeSet<String>();
        props.addAll(p.values().iterator().next().keySet());

        TableWriter tw = new TableWriter(out);

        TableReader tr = new TableReader(CALLS);
        for (Map<String, String> te : tr) {
            String[] locus = te.get("locus").split("[:-]");

            Interval intervalSmall = locus[0].equals("unknown") ? new Interval("*", 0, 0) : new Interval(locus[0], Integer.valueOf(locus[1]), Integer.valueOf(locus[2]));
            Interval interval = locus[0].equals("unknown") ? new Interval("*", 0, 0) : new Interval(locus[0], Integer.valueOf(locus[1]) - 100000, Integer.valueOf(locus[2]) + 100000);

            te.put("closestGene", "NA");
            te.put("withinGene", "NA");
            te.put("withinCDS", "NA");
            te.put("geneDescription", "NA");
            for (String prop : props) { te.put(prop, "NA"); }

            Collection<GFF3Record> genes = new HashSet<GFF3Record>();
            for (GFF3 gff : GFFS) {
                genes.addAll(GFF3.getType("gene", gff.getOverlapping(interval)));
            }

            if (genes.size() > 0) {
                GFF3Record closestGene = genes.iterator().next();
                for (GFF3Record gene : genes) {
                    if (Math.abs(gene.getStart() - intervalSmall.getStart()) < Math.abs(closestGene.getStart() - intervalSmall.getStart())) {
                        closestGene = gene;
                    }
                }

                te.put("closestGene", closestGene.getAttribute("ID"));
                te.put("withinGene", "false");
                te.put("withinCDS", "false");
                te.put("geneDescription", closestGene.getAttribute("description"));

                if (intervalSmall.intersects(closestGene.getInterval())) {
                    te.put("withinGene", "true");
                }

                for (GFF3 gff : GFFS) {
                    for (GFF3Record exon : GFF3.getType("exon", gff.getContained(closestGene))) {
                        if (intervalSmall.intersects(exon.getInterval())) {
                            te.put("withinCDS", "true");
                        }
                    }
                }
            }

            if (p.containsOverlapping(intervalSmall)) {
                Map<String, String> pr = p.getOverlapping(intervalSmall).iterator().next();

                for (String prop : pr.keySet()) {
                    te.put(prop, pr.get(prop));
                }
            }

            tw.addEntry(te);
        }
    }
}
