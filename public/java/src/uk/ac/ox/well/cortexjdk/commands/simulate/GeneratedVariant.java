package uk.ac.ox.well.cortexjdk.commands.simulate;

import org.jetbrains.annotations.NotNull;

public class GeneratedVariant implements Comparable<GeneratedVariant> {
    public String type;
    public int seqIndex;
    public int posIndex;
    public String oldAllele;
    public String newAllele;

    public GeneratedVariant(String type, int seqIndex, int posIndex, String oldAllele, String newAllele) {
        this.type = type;
        this.oldAllele = oldAllele;
        this.newAllele = newAllele;
        this.seqIndex = seqIndex;
        this.posIndex = posIndex;
    }

    public String getType() { return type; }
    public int getSeqIndex() { return seqIndex; }
    public int getPosIndex() { return posIndex; }
    public String getOldAllele() { return oldAllele; }
    public String getNewAllele() { return newAllele; }

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
