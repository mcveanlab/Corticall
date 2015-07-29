package uk.ac.ox.well.indiana.utils.alignment.kmer;

public class GenomePosition {
    public byte chr;
    public int pos;

    public GenomePosition(byte chr, int pos) {
        this.chr = chr;
        this.pos = pos;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        GenomePosition that = (GenomePosition) o;

        if (chr != that.chr) return false;
        return pos == that.pos;

    }

    @Override
    public int hashCode() {
        int result = (int) chr;
        result = 31 * result + pos;
        return result;
    }
}
