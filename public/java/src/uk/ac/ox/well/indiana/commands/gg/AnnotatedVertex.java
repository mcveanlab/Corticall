package uk.ac.ox.well.indiana.commands.gg;

/**
 * Created by kiran on 18/06/2015.
 */
public class AnnotatedVertex {
    private String kmer;
    private boolean isNovel = false;

    public AnnotatedVertex(String kmer) { this.kmer = kmer; }
    public String getKmer() { return kmer; }
    public boolean isNovel() { return isNovel; }
    public void setNovel() { isNovel = true; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        AnnotatedVertex that = (AnnotatedVertex) o;

        return kmer.equals(that.kmer);
    }

    @Override
    public int hashCode() {
        return kmer.hashCode();
    }
}
