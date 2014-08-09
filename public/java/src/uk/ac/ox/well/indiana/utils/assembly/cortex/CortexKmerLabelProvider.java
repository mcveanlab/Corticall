package uk.ac.ox.well.indiana.utils.assembly.cortex;

import org.jgrapht.ext.VertexNameProvider;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

public class CortexKmerLabelProvider implements VertexNameProvider<CortexKmer> {
    @Override
    public String getVertexName(CortexKmer cortexKmer) {
        String fullName = cortexKmer.isFlipped() ? SequenceUtils.reverseComplement(cortexKmer.getKmerAsString()) : cortexKmer.getKmerAsString();

        StringBuilder shortName = new StringBuilder();
        shortName.append(fullName.charAt(0));
        shortName.append(fullName.charAt(1));
        shortName.append(fullName.charAt(2));
        shortName.append("...");
        shortName.append(fullName.charAt(fullName.length() - 3));
        shortName.append(fullName.charAt(fullName.length() - 2));
        shortName.append(fullName.charAt(fullName.length() - 1));

        return shortName.toString();
    }
}
