package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import htsjdk.samtools.util.Interval;
import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;

import java.util.List;

public class GeneratedVariant implements Comparable<GeneratedVariant> {
    public String type;
    public int seqIndex;
    public int posIndex;
    public String oldAllele;
    public String newAllele;
    public List<Pair<Interval, Interval>> loci;

    public GeneratedVariant(String type, int seqIndex, int posIndex, String oldAllele, String newAllele) {
        this.type = type;
        this.oldAllele = oldAllele;
        this.newAllele = newAllele;
        this.seqIndex = seqIndex;
        this.posIndex = posIndex;
    }

    public GeneratedVariant(String type, int seqIndex, int posIndex, String oldAllele, String newAllele, List<Pair<Interval, Interval>> loci) {
        this.type = type;
        this.oldAllele = oldAllele;
        this.newAllele = newAllele;
        this.seqIndex = seqIndex;
        this.posIndex = posIndex;
        this.loci = loci;
    }

    public String getType() { return type; }
    public int getSeqIndex() { return seqIndex; }
    public int getPosIndex() { return posIndex; }
    public String getOldAllele() { return oldAllele; }
    public String getNewAllele() { return newAllele; }
    public List<Pair<Interval, Interval>> getLoci() { return loci; }

    @Override
    public int compareTo(@NotNull GeneratedVariant o) {
        if (getSeqIndex() != o.getSeqIndex()) {
            return getSeqIndex() < o.getSeqIndex() ? -1 : 1;
        }

        if (getPosIndex() != o.getPosIndex()) {
            return getPosIndex() < o.getPosIndex() ? -1 : 1;
        }

        return 0;
    }

    @Override
    public String toString() {
        return "GeneratedVariant{" +
                "type='" + type + '\'' +
                ", seqIndex=" + seqIndex +
                ", posIndex=" + posIndex +
                ", oldAllele='" + oldAllele + '\'' +
                ", newAllele='" + newAllele + '\'' +
                '}';
    }
}
