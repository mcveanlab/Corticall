package uk.ac.ox.well.cortexjdk.utils.io.graph;

import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexHeader;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinksRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexBinaryKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;

import java.io.File;

/**
 * Created by kiran on 14/09/2017.
 */
public interface ConnectivityAnnotations {
    File getFile();

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
        } else if (key instanceof CanonicalKmer) {
            return new CortexBinaryKmer(((CanonicalKmer) key).getKmerAsBytes());
        } else if (key instanceof String) {
            return new CortexBinaryKmer(((String) key).getBytes());
        } else if (key instanceof CortexByteKmer) {
            return new CortexBinaryKmer(((CortexByteKmer) key).getKmer());
        }

        throw new CortexJDKException("Could not convert object to internal representation (" + key + ")");
    }
}
