package uk.ac.ox.well.cortexjdk.utils.caller;

import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.VariantContext;
import org.jgrapht.GraphPath;
import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;

import java.util.Set;

public class Bubble {
    private String flank5p;
    private String refAllele;
    private String altAllele;
    private String flank3p;

    private Set<CortexKmer> novelKmers;
    private Interval locus;

    public Bubble(GraphPath<CortexVertex, CortexEdge> pRef, GraphPath<CortexVertex, CortexEdge> pAlt, KmerLookup kl, Set<CortexKmer> novelKmers) {
        String[] pieces = pathsToAlleles(pRef, pAlt);

        this.flank5p = pieces[0];
        this.refAllele = pieces[1];
        this.altAllele = pieces[2];
        this.flank3p = pieces[3];
        this.novelKmers = novelKmers;

        Interval li5p = null, li3p = null;

        for (int i = 0; i <= flank5p.length() - kl.getKmerSize(); i++) {
            String sk = flank5p.substring(i, i + kl.getKmerSize());
            Set<Interval> lis = kl.findKmer(sk);

            if (lis.size() == 1) {
                Interval li = lis.iterator().next();
                int offset = flank5p.length() - i;

                if (li.isNegativeStrand()) {
                    offset += refAllele.length();
                    offset *= -1;
                }

                String contig = li.getContig();
                int pos = li.getStart() + offset;

                li5p = new Interval(contig, pos, pos);
                break;
            }
        }

        for (int i = 0; i <= flank3p.length() - kl.getKmerSize(); i++) {
            String sk = flank3p.substring(i, i + kl.getKmerSize());
            Set<Interval> lis = kl.findKmer(sk);

            if (lis.size() == 1) {
                Interval li = lis.iterator().next();
                int offset = i;

                if (li.isPositiveStrand()) {
                    offset += refAllele.length();
                    offset *= -1;
                }

                String contig = li.getContig();
                int pos = li.getStart() + offset;

                li3p = new Interval(contig, pos, pos);
                break;
            }
        }

        if (li5p != null) {
            locus = li5p;
        } else if (li3p != null) {
            locus = li3p;
        }
    }

    public VariantContext.Type getType() {
        if (refAllele.length() == altAllele.length() && !refAllele.equals(altAllele)) {
            if (refAllele.length() == 1) {
                return VariantContext.Type.SNP;
            }

            return VariantContext.Type.MNP;
        } else if (refAllele.length() != altAllele.length()) {
            if (refAllele.length() == 0 || altAllele.length() == 0) {
                return VariantContext.Type.INDEL;
            }

            return VariantContext.Type.MNP;
        }

        return VariantContext.Type.NO_VARIATION;
    }

    public String getFlank5p() { return flank5p; }

    public String getFlank3p() { return flank3p; }

    public String getRefAllele() {
        return refAllele;
    }

    public String getAltAllele() {
        return altAllele;
    }

    public String getRefHaplotype() { return refAllele; }

    public String getAltHaplotype() {
        return altAllele;
    }

    public Set<CortexKmer> getNovelKmers() { return novelKmers; }

    public Interval getLocus() { return locus; }

    @Override
    public String toString() {
        return "Bubble{" +
                "refAllele='" + refAllele + '\'' +
                ", altAllele='" + altAllele + '\'' +
                ", locus='" + locus + '\'' +
                '}';
    }

    private static String[] pathsToAlleles(GraphPath<CortexVertex, CortexEdge> p0, GraphPath<CortexVertex, CortexEdge> p1) {
        String s0 = pathToString(p0);
        String s1 = pathToString(p1);

        int s0start = 0, s0end = s0.length();
        int s1start = 0, s1end = s1.length();

        for (int i = 0, j = 0; i < s0.length() && j < s1.length(); i++, j++) {
            if (s0.charAt(i) != s1.charAt(j)) {
                s0start = i;
                s1start = j;
                break;
            }
        }

        for (int i = s0.length() - 1, j = s1.length() - 1; i >= 0 && j >= 0; i--, j--) {
            if (s0.charAt(i) != s1.charAt(j) || i == s0start - 1 || j == s1start - 1) {
                s0end = i + 1;
                s1end = j + 1;
                break;
            }
        }

        String[] pieces = new String[4];
        pieces[0] = s0.substring(0, s0start);
        pieces[1] = s0.substring(s0start, s0end);
        pieces[2] = s1.substring(s1start, s1end);
        pieces[3] = s0.substring(s0end, s0.length());

        return pieces;
    }

    private static String pathToString(GraphPath<CortexVertex, CortexEdge> pk) {
        StringBuilder sbk = new StringBuilder();
        if (pk != null) {
            for (CortexVertex v : pk.getVertexList()) {
                if (sbk.length() == 0) {
                    sbk.append(v.getSk());
                } else {
                    sbk.append(v.getSk().substring(v.getSk().length() - 1, v.getSk().length()));
                }
            }
        }
        return sbk.toString();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Bubble bubble = (Bubble) o;

        if (refAllele != null ? !refAllele.equals(bubble.refAllele) : bubble.refAllele != null) return false;
        if (altAllele != null ? !altAllele.equals(bubble.altAllele) : bubble.altAllele != null) return false;
        if (novelKmers != null ? !novelKmers.equals(bubble.novelKmers) : bubble.novelKmers != null) return false;
        return locus != null ? locus.equals(bubble.locus) : bubble.locus == null;

    }

    @Override
    public int hashCode() {
        int result = refAllele != null ? refAllele.hashCode() : 0;
        result = 31 * result + (altAllele != null ? altAllele.hashCode() : 0);
        result = 31 * result + (novelKmers != null ? novelKmers.hashCode() : 0);
        result = 31 * result + (locus != null ? locus.hashCode() : 0);
        return result;
    }
}
