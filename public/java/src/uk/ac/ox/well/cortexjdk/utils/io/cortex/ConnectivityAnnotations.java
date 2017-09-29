package uk.ac.ox.well.cortexjdk.utils.io.cortex;

import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexBinaryKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexHeader;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinksRecord;

/**
 * Created by kiran on 14/09/2017.
 */
public interface ConnectivityAnnotations {
    int size();

    boolean isEmpty();

    boolean containsKey(Object key);

    CortexLinksRecord get(Object key);

    CortexHeader getHeader();

    default String getSource() { return "unknown"; }

    default CortexBinaryKmer convert(Object key) {
        if (key instanceof CortexBinaryKmer) {
            return (CortexBinaryKmer) key;
        } else if (key instanceof byte[]) {
            return new CortexBinaryKmer((byte[]) key);
        } else if (key instanceof CortexKmer) {
            return new CortexBinaryKmer(((CortexKmer) key).getKmerAsBytes());
        } else if (key instanceof String) {
            return new CortexBinaryKmer(((String) key).getBytes());
        } else if (key instanceof CortexByteKmer) {
            return new CortexBinaryKmer(((CortexByteKmer) key).getKmer());
        }

        throw new CortexJDKException("Could not convert object to internal representation (" + key + ")");
    }
}
