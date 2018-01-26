package uk.ac.ox.well.cortexjdk.utils.traversal;

import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;

import java.util.HashSet;
import java.util.Set;

public class CortexVertexFactory {
    private CortexByteKmer bk;
    private CortexRecord cr;
    private Interval locus = null;
    private Set<String> kmerSources = new HashSet<>();
    private int copyIndex = 0;

    public CortexVertexFactory vertex(CortexVertex v) {
        bases(v.getKmerAsByteKmer());
        record(v.getCortexRecord());
        locus(v.getLocus());
        sources(v.getSources());
        copyIndex(v.getCopyIndex());

        return this;
    }

    public CortexVertexFactory bases(String sk) { this.bk = new CortexByteKmer(sk); return this; }
    public CortexVertexFactory bases(CortexByteKmer bk) { this.bk = bk; return this; }
    public CortexVertexFactory bases(byte[] bk) { this.bk = new CortexByteKmer(bk); return this; }

    public CortexVertexFactory record(CortexRecord cr) { this.cr = cr; return this; }

    public CortexVertexFactory locus(Interval locus) { this.locus = locus; return this; }

    public CortexVertexFactory sources(Set<String> sources) { this.kmerSources = sources; return this; }
    public CortexVertexFactory source(String source) { this.kmerSources.add(source); return this; }

    public CortexVertexFactory copyIndex(int copyIndex) { this.copyIndex = copyIndex; return this; }

    public CortexVertex make() {
        return new CortexVertex(bk, cr, locus, kmerSources, copyIndex);
    }
}
