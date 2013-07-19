package uk.ac.ox.well.indiana.utils.assembly;

import org.jgrapht.ext.VertexNameProvider;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

public class CortexKmerIDProvider implements VertexNameProvider<CortexKmer> {
    @Override
    public String getVertexName(CortexKmer cortexKmer) {
        return cortexKmer.isFlipped() ? SequenceUtils.reverseComplement(cortexKmer.getKmerAsString()) : cortexKmer.getKmerAsString();
    }
}
