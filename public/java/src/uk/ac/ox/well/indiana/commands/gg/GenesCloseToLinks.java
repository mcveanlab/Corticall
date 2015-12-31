package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class GenesCloseToLinks extends Module {
    @Argument(fullName="links", shortName="l", doc="Links")
    public File LINKS;

    @Argument(fullName="gffs", shortName="g", doc="Gff")
    public ArrayList<GFF3> GFFS;

    @Argument(fullName="proximity", shortName="p", doc="Gene proximity")
    public Integer PROXIMITY = 50;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Set<GFF3Record> records = new HashSet<GFF3Record>();

        LineReader lr = new LineReader(LINKS);

        Map<String, Boolean> overlaps = new HashMap<String, Boolean>();

        while (lr.hasNext()) {
            String[] p = lr.getNextRecord().split("\\s+");

            String chr = p[1];
            int start = Integer.valueOf(p[2]);
            int stop = Integer.valueOf(p[3]);

            Interval interval = new Interval(chr, start - 100000, stop + 100000);

            Collection<GFF3Record> genes = new HashSet<GFF3Record>();
            for (GFF3 gff : GFFS) {
                genes.addAll(GFF3.getType("gene", gff.getOverlapping(interval)));
            }

            if (genes.size() > 0) {
                GFF3Record closestGene = genes.iterator().next();
                for (GFF3Record gene : genes) {
                    if (Math.abs(gene.getStart() - start) < Math.abs(closestGene.getStart() - start)) {
                        closestGene = gene;
                    }
                }

                if (start >= closestGene.getStart() && stop <= closestGene.getEnd()) {
                    overlaps.put(closestGene.getAttribute("ID"), true);
                } else {
                    overlaps.put(closestGene.getAttribute("ID"), false);
                }

                log.info("  {} {}:{}-{} {} {}", closestGene.getAttribute("ID"), closestGene.getSeqid(), closestGene.getStart(), closestGene.getEnd(), start, stop);

                records.add(closestGene);
            }
        }

        for (GFF3Record g : records) {
            String style = "color=" + (overlaps.get(g.getAttribute("ID")) ? "black" : "grey");
            String name = g.getAttribute("ID");
            out.println(Joiner.on(" ").join(g.getSeqid(), g.getStart(), g.getEnd(), name, style));
        }
    }
}
