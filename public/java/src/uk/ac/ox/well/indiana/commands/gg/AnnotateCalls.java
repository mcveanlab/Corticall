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
            List<String> closestGenes = new ArrayList<String>();
            List<String> withinGenes = new ArrayList<String>();
            List<String> withinCDS = new ArrayList<String>();
            List<String> geneDescription = new ArrayList<String>();

            for (String prop : props) {
                te.put(prop, "NA");
            }

            if (!te.get("locus").contains("unknown")) {
                String[] loci = te.get("locus").split(";");

                for (String locus : loci) {
                    String[] pieces = locus.split("[:-]");

                    Interval intervalSmall = new Interval(pieces[0], Integer.valueOf(pieces[1]),          Integer.valueOf(pieces[2]));
                    Interval interval      = new Interval(pieces[0], Integer.valueOf(pieces[1]) - 100000, Integer.valueOf(pieces[2]) + 100000);

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

                        String cg = closestGene.getAttribute("ID");
                        boolean inGene = false;
                        boolean inCDS = false;
                        String desc = closestGene.getAttribute("description");

                        if (intervalSmall.intersects(closestGene.getInterval())) {
                            inGene = true;
                        }

                        for (GFF3 gff : GFFS) {
                            for (GFF3Record exon : GFF3.getType("exon", gff.getContained(closestGene))) {
                                if (intervalSmall.intersects(exon.getInterval())) {
                                    inCDS = true;
                                }
                            }
                        }

                        closestGenes.add(cg);
                        withinGenes.add(inGene ? "true" : "false");
                        withinCDS.add(inCDS ? "true" : "false");
                        geneDescription.add(desc);
                    }

                    if (p.containsOverlapping(intervalSmall)) {
                        Map<String, String> pr = p.getOverlapping(intervalSmall).iterator().next();

                        for (String prop : pr.keySet()) {
                            te.put(prop, pr.get(prop));
                        }
                    }
                }

                te.put("closestGene", Joiner.on(";").join(closestGenes));
                te.put("withinGene", Joiner.on(";").join(withinGenes));
                te.put("withinCDS", Joiner.on(";").join(withinCDS));
                te.put("geneDescription", Joiner.on(";").join(geneDescription));
            } else {
                te.put("closestGene", "unknown");
                te.put("withinGene", "NA");
                te.put("withinCDS", "NA");
                te.put("geneDescription", "NA");
            }

            tw.addEntry(te);
        }
    }
}
