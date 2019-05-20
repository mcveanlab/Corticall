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
    public String seedLeft;
    public String seedRight;
    public String parent;
    public Interval start;
    public Interval stop;

    public GeneratedVariant(String type, int seqIndex, int posIndex, String oldAllele, String newAllele) {
        this.type = type;
        this.oldAllele = oldAllele;
        this.newAllele = newAllele;
        this.seqIndex = seqIndex;
        this.posIndex = posIndex;
    }

    public GeneratedVariant(String type, int seqIndex, int posIndex, String oldAllele, String newAllele, String seedLeft, String seedRight) {
        this.type = type;
        this.oldAllele = oldAllele;
        this.newAllele = newAllele;
        this.seqIndex = seqIndex;
        this.posIndex = posIndex;
        this.seedLeft = seedLeft;
        this.seedRight = seedRight;
    }

    public GeneratedVariant(String type, int seqIndex, int posIndex, String oldAllele, String newAllele, String seedLeft, String seedRight, String parent, Interval start, Interval stop) {
        this.type = type;
        this.oldAllele = oldAllele;
        this.newAllele = newAllele;
        this.seqIndex = seqIndex;
        this.posIndex = posIndex;
        this.seedLeft = seedLeft;
        this.seedRight = seedRight;
        this.parent = parent;
        this.start = start;
        this.stop = stop;
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

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        GeneratedVariant that = (GeneratedVariant) o;

        if (seqIndex != that.seqIndex) return false;
        if (posIndex != that.posIndex) return false;
        if (type != null ? !type.equals(that.type) : that.type != null) return false;
        if (oldAllele != null ? !oldAllele.equals(that.oldAllele) : that.oldAllele != null) return false;
        if (newAllele != null ? !newAllele.equals(that.newAllele) : that.newAllele != null) return false;
        if (loci != null ? !loci.equals(that.loci) : that.loci != null) return false;
        if (seedLeft != null ? !seedLeft.equals(that.seedLeft) : that.seedLeft != null) return false;
        if (seedRight != null ? !seedRight.equals(that.seedRight) : that.seedRight != null) return false;
        return parent != null ? parent.equals(that.parent) : that.parent == null;
    }

    @Override
    public int hashCode() {
        int result = type != null ? type.hashCode() : 0;
        result = 31 * result + seqIndex;
        result = 31 * result + posIndex;
        result = 31 * result + (oldAllele != null ? oldAllele.hashCode() : 0);
        result = 31 * result + (newAllele != null ? newAllele.hashCode() : 0);
        result = 31 * result + (loci != null ? loci.hashCode() : 0);
        result = 31 * result + (seedLeft != null ? seedLeft.hashCode() : 0);
        result = 31 * result + (seedRight != null ? seedRight.hashCode() : 0);
        result = 31 * result + (parent != null ? parent.hashCode() : 0);
        return result;
    }
}
