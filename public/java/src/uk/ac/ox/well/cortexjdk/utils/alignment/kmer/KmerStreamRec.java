package uk.ac.ox.well.cortexjdk.utils.alignment.kmer;

import uk.ac.ox.well.cortexjdk.utils.kmer.CortexBinaryKmer;

/**
 * Created by kiran on 20/09/2017.
 */
public class KmerStreamRec {
    public CortexBinaryKmer cbk;
    public int contigIndex;
    public int offset;

    public KmerStreamRec(long[] l, int ci, int o) {
        cbk = new CortexBinaryKmer(l);
        contigIndex = ci;
        offset = o;
    }

    public int kmerCompare(KmerStreamRec that) {
        return cbk.compareTo(that.cbk);
    }
}

