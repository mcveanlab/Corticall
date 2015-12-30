package uk.ac.ox.well.indiana.commands.gg;

import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

import java.io.File;
import java.io.PrintStream;
import java.util.Map;

public class AnnotateCalls extends Module {
    @Argument(fullName="links", shortName="l", doc="Circos links")
    public File LINKS;

    @Argument(fullName="gff1", shortName="g1", doc="GFF1")
    public GFF3 GFF1;

    @Argument(fullName="gff2", shortName="g2", doc="GFF2")
    public GFF3 GFF2;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        //PG0030-C.ERR019046_stretch945 PfDd2_13_TT 2060264 2062997 radius2=0.6r # unknown_incomplete
        TableReader tr = new TableReader(LINKS, "id", "chr", "start", "stop", "radius", "comment", "type");

        for (Map<String, String> te : tr) {
            Interval i = new Interval(te.get("chr"), Integer.valueOf(te.get("start")), Integer.valueOf(te.get("stop")));

            log.info("i 1 {} {}", i, GFF1.getOverlapping(i));
            log.info("i 2 {} {}", i, GFF2.getOverlapping(i));
        }
    }
}
