package uk.ac.ox.well.indiana.utils.assembly.cortex;

import org.jgrapht.ext.ComponentAttributeProvider;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;

import java.util.HashMap;
import java.util.Map;

public class CortexKmerAttributeProvider implements ComponentAttributeProvider<CortexKmer> {
    private CortexKmer panelKmer;

    public CortexKmerAttributeProvider(CortexKmer panelKmer) {
        this.panelKmer = panelKmer;
    }

    @Override
    public Map<String, String> getComponentAttributes(CortexKmer cortexKmer) {
        Map<String, String> attrs = new HashMap<String, String>();

        if (panelKmer.equals(cortexKmer)) {
            attrs.put("color", "red");
        }

        return attrs;
    }
}
