package uk.ac.ox.well.cortexjdk.utils.traversal;

import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;

import java.util.Set;

/**
 * Created by kiran on 10/05/2017.
 */
public class CortexVertex {
    private CortexByteKmer sk;
    private CortexRecord cr;
    private Interval locus;
    private Set<String> kmerSources;
    private int copyIndex = 0;

    public CortexVertex(CortexByteKmer sk, CortexRecord cr) {
        this.sk = sk;
        this.cr = cr;
    }

    public CortexVertex(CortexByteKmer sk, CortexRecord cr, Interval locus) {
        this.sk = sk;
        this.cr = cr;
        this.locus = locus;
    }

    public CortexVertex(CortexByteKmer sk, CortexRecord cr, Set<String> kmerSources) {
        this.sk = sk;
        this.cr = cr;
        this.kmerSources = kmerSources;
    }

    public CortexVertex(CortexByteKmer sk, CortexRecord cr, Interval locus, Set<String> kmerSources) {
        this.sk = sk;
        this.cr = cr;
        this.locus = locus;
        this.kmerSources = kmerSources;
    }

    public CortexVertex(CortexByteKmer sk, CortexRecord cr, Interval locus, Set<String> kmerSources, int copyIndex) {
        this.sk = sk;
        this.cr = cr;
        this.locus = locus;
        this.kmerSources = kmerSources;
        this.copyIndex = copyIndex;
    }

    public String getKmerAsString() { return new String(sk.getKmer()); }

    public CortexByteKmer getKmerAsByteKmer() { return sk; }

    public CortexRecord getCortexRecord() { return cr; }

    public CanonicalKmer getCanonicalKmer() { return cr.getCanonicalKmer(); }

    public int getCopyIndex() { return copyIndex; }

    public Interval getLocus() { return locus; }

    public Set<String> getSources() { return kmerSources; }

    public void setLocus(Interval locus) {
        this.locus = locus;
    }

    public void setKmerSources(Set<String> kmerSources) {
        this.kmerSources = kmerSources;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        CortexVertex that = (CortexVertex) o;

        if (copyIndex != that.copyIndex) return false;
        if (!sk.equals(that.sk)) return false;
        if (!cr.equals(that.cr)) return false;
        if (locus != null ? !locus.equals(that.locus) : that.locus != null) return false;
        return kmerSources != null ? kmerSources.equals(that.kmerSources) : that.kmerSources == null;
    }

    @Override
    public int hashCode() {
        int result = sk.hashCode();
        result = 31 * result + cr.hashCode();
        result = 31 * result + (locus != null ? locus.hashCode() : 0);
        result = 31 * result + (kmerSources != null ? kmerSources.hashCode() : 0);
        result = 31 * result + copyIndex;
        return result;
    }

    @Override
    public String toString() {
        return "CortexVertex{" +
                "sk=" + sk +
                ", cr=" + cr +
                ", locus=" + locus +
                ", kmerSources=" + kmerSources +
                ", copyIndex=" + copyIndex +
                '}';
    }
}
