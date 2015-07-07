package uk.ac.ox.well.indiana.commands.gg;

/**
 * Created by kiran on 18/06/2015.
 */
public class AnnotatedVertex {
    private String kmer;
    private boolean isNovel = false;

    public AnnotatedVertex(String kmer) {
        this.kmer = kmer;
        this.isNovel = false;
    }

    public AnnotatedVertex(String kmer, boolean novelty) {
        this.kmer = kmer;
        this.isNovel = novelty;
    }

    public String getKmer() { return kmer; }
    public boolean isNovel() { return isNovel; }
    public void setNovel() { isNovel = true; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        AnnotatedVertex that = (AnnotatedVertex) o;

        if (isNovel != that.isNovel) return false;
        return !(kmer != null ? !kmer.equals(that.kmer) : that.kmer != null);
    }

    @Override
    public int hashCode() {
        int result = kmer != null ? kmer.hashCode() : 0;
        result = 31 * result + (isNovel ? 1 : 0);
        return result;
    }

    @Override
    public String toString() {
        return "AnnotatedVertex{" +
                "kmer='" + kmer + '\'' +
                ", isNovel=" + isNovel +
                '}';
    }
}
