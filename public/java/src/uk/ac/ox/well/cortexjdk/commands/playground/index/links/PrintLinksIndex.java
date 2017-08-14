package uk.ac.ox.well.cortexjdk.commands.playground.index.links;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinks;

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
