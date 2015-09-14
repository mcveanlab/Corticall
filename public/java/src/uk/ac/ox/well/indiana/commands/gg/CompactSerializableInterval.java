package uk.ac.ox.well.indiana.commands.gg;

import htsjdk.samtools.util.Interval;

import java.io.Serializable;

public class CompactSerializableInterval implements Serializable {
    private int chrIndex;
    private int start;

    public CompactSerializableInterval() {}

    public CompactSerializableInterval(int chrIndex, int start) {
        this.chrIndex = chrIndex;
        this.start = start;
    }

    public void setChrIndex(int chrIndex) { this.chrIndex = chrIndex; }
    public void setStart(int start) { this.start = start; }

    public int getChrIndex() { return chrIndex; }
    public int getStart() { return start; }

    @Override
    public String toString() {
        return "CSI{" +
                "c=" + chrIndex +
                ", s=" + start +
                '}';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        CompactSerializableInterval that = (CompactSerializableInterval) o;

        if (chrIndex != that.chrIndex) return false;
        return start == that.start;

    }

    @Override
    public int hashCode() {
        int result = chrIndex;
        result = 31 * result + start;
        return result;
    }
}
