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

    @Argument(fullName="gff", shortName="g", doc="Gff")
    public GFF3 GFF;

    @Argument(fullName="proximity", shortName="p", doc="Gene proximity")
    public Integer PROXIMITY = 500;

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

            for (int w = 0; w < 1000000; w++) {
                Interval interval = new Interval(chr, start - w, stop + w);

                Collection<GFF3Record> genes = GFF3.getType("gene", GFF.getOverlapping(interval));

                if (genes.size() > 0) {
                    records.addAll(genes);

                    for (GFF3Record gr : genes) {
                        if (w < PROXIMITY) {
                            overlaps.put(gr.getAttribute("ID"), true);
                        } else {
                            overlaps.put(gr.getAttribute("ID"), false);
                        }

                        log.info("  {}:{}-{} {}", gr.getSeqid(), gr.getStart(), gr.getEnd(), w);
                    }

                    break;
                }
            }
        }

        for (GFF3Record g : records) {
            String style = "color=" + (overlaps.get(g.getAttribute("ID")) ? "black" : "red");
            String name = g.getAttribute("ID");
            out.println(Joiner.on(" ").join(g.getSeqid(), g.getStart(), g.getEnd(), name, style));
        }
    }
}
