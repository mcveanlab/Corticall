package uk.ac.ox.well.indiana.utils.caller;

import com.google.common.base.Joiner;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.GraphPath;
import uk.ac.ox.well.indiana.utils.alignment.sw.SmithWaterman;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;

import java.util.Set;

public class Bubble {
    private int contigPos;
    private Set<CortexKmer> usedNovelKmers;

    private String flank5p;
    private String refAllele;
    private String altAllele;
    private String flank3p;

    private int refColor;
    private int altColor;

    public Bubble(String flank5p, String refAllele, String altAllele, String flank3p) {
        this.flank5p = flank5p;
        this.refAllele = refAllele;
        this.altAllele = altAllele;
        this.flank3p = flank3p;
    }

    public Bubble(GraphPath<CortexVertex, CortexEdge> pRef, GraphPath<CortexVertex, CortexEdge> pAlt) {
        String[] pieces = pathsToAlleles(pRef, pAlt);

        this.flank5p = pieces[0];
        this.refAllele = pieces[1];
        this.altAllele = pieces[2];
        this.flank3p = pieces[3];
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Bubble bubble = (Bubble) o;

        CortexKmer thisRefHaplotypeFw = new CortexKmer(getRefHaplotype());
        CortexKmer thisAltHaplotypeFw = new CortexKmer(getAltHaplotype());

        CortexKmer thatRefHaplotypeFw = new CortexKmer(bubble.getRefHaplotype());
        CortexKmer thatAltHaplotypeFw = new CortexKmer(bubble.getAltHaplotype());

        return thisRefHaplotypeFw.equals(thatRefHaplotypeFw) && thisAltHaplotypeFw.equals(thatAltHaplotypeFw);
    }

    public String getRefAllele() {
        return refAllele;
    }

    public String getAltAllele() {
        return altAllele;
    }

    public String getRefHaplotype() {
        //return flank5p + refAllele + flank3p;
        return refAllele;
    }

    public String getAltHaplotype() {
        return altAllele;
    }

    @Override
    public int hashCode() {
        int result = refAllele != null ? (new CortexKmer(refAllele)).hashCode() : 0;
        result = 31 * result + (altAllele != null ? (new CortexKmer(altAllele)).hashCode() : 0);
        return result;
    }

    @Override
    public String toString() {
        return "Bubble{" +
                "refAllele='" + refAllele + '\'' +
                ", altAllele='" + altAllele + '\'' +
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
}
