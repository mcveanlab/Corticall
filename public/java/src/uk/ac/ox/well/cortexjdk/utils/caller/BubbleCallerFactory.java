package uk.ac.ox.well.cortexjdk.utils.caller;

import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinks;

import java.util.Arrays;
import java.util.Collection;
import java.util.Map;

/**
 * Created by kiran on 30/08/2017.
 */
public class BubbleCallerFactory {
    private BubbleCallerConfiguration configuration = new BubbleCallerConfiguration();

    public BubbleCallerFactory alternateColor(int color) { configuration.setAlternateColor(color); return this; }

    public BubbleCallerFactory referenceColors(int... colors) { Arrays.stream(colors).forEach(c -> configuration.getReferenceColors().add(c)); return this; }
    public BubbleCallerFactory referenceColors(Collection<Integer> colors) { configuration.getReferenceColors().addAll(colors); return this; }

    public BubbleCallerFactory references(Map<String, KmerLookup> references) { configuration.setReferences(references); return this; }

    public BubbleCallerFactory graph(DeBruijnGraph graph) { configuration.setGraph(graph); return this; }
    public BubbleCallerFactory rois(DeBruijnGraph rois) { configuration.setRois(rois); return this; }
    public BubbleCallerFactory links(CortexLinks... links) { Arrays.stream(links).forEach(l -> configuration.getLinks().add(l)); return this; }
    public BubbleCallerFactory links(Collection<CortexLinks> links) { configuration.getLinks().addAll(links); return this; }

    public BubbleCaller make() {
        return new BubbleCaller(configuration);
    }
}
