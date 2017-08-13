package uk.ac.ox.well.indiana.commands.playground.index.links;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksRecord;

import java.io.PrintStream;

/**
 * Created by kiran on 13/08/2017.
 */
public class PrintLinksIndex extends Module {
    @Argument(fullName="links", shortName="l", doc="Links")
    public CortexLinks LINKS;

    @Override
    public void execute() {
        for (CortexKmer ck : LINKS.keySet()) {
            log.info("{} {}", ck, LINKS.get(ck));
        }
    }
}
