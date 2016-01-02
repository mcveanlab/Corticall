package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
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
    @Argument(fullName="links", shortName="l", doc="Circos links")
    public File LINKS;

    @Argument(fullName="gffs", shortName="g", doc="GFFs")
    public ArrayList<GFF3> GFFS;

    @Argument(fullName="properties", shortName="p", doc="Sequence properties")
    public ArrayList<File> PROPS;

    @Argument(fullName="calls", shortName="d", doc="Calls")
    public File CALLS;

    @Argument(fullName="cross", shortName="c", doc="Cross")
    public String CROSS;

    @Output
    public PrintStream out;

    private Map<String, String> parseLine(String line) {
        String[] pieces = line.split("\\s+");

        Map<String, String> te = new HashMap<String, String>();
        te.put("s_id", pieces[0]);
        te.put("chr", pieces[1]);
        te.put("start", pieces[2]);
        te.put("stop", pieces[3]);
        te.put("radius", pieces[4]);
        te.put("comment", pieces[5]);
        te.put("type_completeness", pieces[6]);

        return te;
    }

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

    private Map<String, Map<String, String>> loadCalls() {
        Map<String, Map<String, String>> calls = new HashMap<String, Map<String, String>>();

        TableReader tr = new TableReader(CALLS);

        for (Map<String, String> te : tr) {
            String key = String.format("%s%s%s", te.get("cross"), te.get("sample"), te.get("id"));

            calls.put(key, te);
        }

        return calls;
    }

    @Override
    public void execute() {
        Set<Map<String, String>> entries = new LinkedHashSet<Map<String, String>>();

        IntervalTreeMap<Map<String, String>> p = loadProperties();
        Map<String, Map<String, String>> c = loadCalls();

        Set<String> props = new TreeSet<String>();
        props.addAll(p.values().iterator().next().keySet());

        LineReader lr = new LineReader(LINKS);
        while (lr.hasNext()) {
            Map<String, String> te = parseLine(lr.getNextRecord());

            String[] sampleAndId = te.get("s_id").split("_stretch");
            //String[] sampleAndId = te.get("s_id").split("_", 1);

            String[] typeAndCompleteness = te.get("type_completeness").replaceAll("STR_", "STR").split("_");

            Interval intervalSmall = new Interval(te.get("chr"), Integer.valueOf(te.get("start")), Integer.valueOf(te.get("stop")));
            Interval interval = new Interval(intervalSmall.getSequence(), intervalSmall.getStart() - 100000, intervalSmall.getEnd() + 100000);

            Map<String, String> entry = new LinkedHashMap<String, String>();
            entry.put("cross", CROSS);
            entry.put("sample", sampleAndId[0]);
            entry.put("id", "stretch" + sampleAndId[1].split("\\.")[0]);
            entry.put("subid", "stretch" + sampleAndId[1]);
            entry.put("type", typeAndCompleteness[0].replace("STR", "STR_"));
            entry.put("completeness", typeAndCompleteness[1]);
            entry.put("locus", String.format("%s:%d-%d", intervalSmall.getSequence(), intervalSmall.getStart(), intervalSmall.getEnd()));
            entry.put("closestGene", "NA");
            entry.put("withinGene", "NA");
            entry.put("withinCDS", "NA");
            entry.put("geneDescription", "NA");
            for (String prop : props) {
                entry.put(prop, "NA");
            }

            entry.put("type", "NA");
            entry.put("completeness", "NA");
            entry.put("parentalAllele", "NA");
            entry.put("childAllele", "NA");
            entry.put("novelKmer", "NA");

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

                entry.put("closestGene", closestGene.getAttribute("ID"));
                entry.put("withinGene", "false");
                entry.put("withinCDS", "false");
                entry.put("geneDescription", closestGene.getAttribute("description"));

                if (intervalSmall.intersects(closestGene.getInterval())) {
                    entry.put("withinGene", "true");
                }

                for (GFF3 gff : GFFS) {
                    for (GFF3Record exon : GFF3.getType("exon", gff.getContained(closestGene))) {
                        if (intervalSmall.intersects(exon.getInterval())) {
                            entry.put("withinCDS", "true");
                        }
                    }
                }
            }

            if (p.containsOverlapping(intervalSmall)) {
                Map<String, String> pr = p.getOverlapping(intervalSmall).iterator().next();

                for (String prop : pr.keySet()) {
                    entry.put(prop, pr.get(prop));
                }
            }

            String key = String.format("%s%s%s", entry.get("cross"), entry.get("sample"), entry.get("id"));

            if (c.containsKey(key)) {
                for (String prop : c.get(key).keySet()) {
                    if (!prop.equals("cross") && !prop.equals("sample") && !prop.equals("id")) {
                        entry.put(prop, c.get(key).get(prop));
                    }
                }
            }

            entries.add(entry);
        }

        TableWriter tw = new TableWriter(out);

        for (Map<String, String> entry : entries) {
            tw.addEntry(entry);
        }

    }
}
