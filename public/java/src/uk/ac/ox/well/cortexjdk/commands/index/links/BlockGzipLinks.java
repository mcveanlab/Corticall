package uk.ac.ox.well.cortexjdk.commands.index.links;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Description;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinksIterable;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinksRecord;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

@Description(text="Block-compresses a .ctp.gz file")
public class BlockGzipLinks extends Module {
    @Argument(fullName="links", shortName="l", doc="Links")
    public CortexLinksIterable LINKS;

    @Override
    public void execute() {
        Path bgzipPath = Paths.get(LINKS.getCortexLinksFile().getAbsolutePath().replace(".ctp.gz", ".ctp.bgz"));

        try {
            BlockCompressedOutputStream bc = new BlockCompressedOutputStream(bgzipPath.toFile().getAbsolutePath(), 9);

            ProgressMeter pm = new ProgressMeterFactory()
                    .header("Processing links")
                    .message("links")
                    .maxRecord(LINKS.getNumKmersWithLinks())
                    .make(log);

            bc.write(LINKS.getHeader().getBytes());
            bc.write("\n".getBytes());
            bc.write(LINKS.getComments().getBytes());
            bc.write("\n".getBytes());

            for (CortexLinksRecord clr : LINKS) {
                bc.write(clr.toString().getBytes());

                pm.update();
            }

            bc.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
