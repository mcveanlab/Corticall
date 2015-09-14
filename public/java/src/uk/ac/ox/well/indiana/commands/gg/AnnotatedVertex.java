package uk.ac.ox.well.indiana.commands.gg;


import htsjdk.samtools.util.Interval;

import java.util.Set;
import java.util.TreeSet;

public class AnnotatedVertex {
    private String kmer;
    private boolean isNovel = false;
    private Set<String> flags = new TreeSet<String>();
    private Set<Interval> ml = new TreeSet<Interval>();
    private Set<Interval> pl = new TreeSet<Interval>();

    public AnnotatedVertex(String kmer) {
        this.kmer = kmer;
        this.isNovel = false;
    }

    public AnnotatedVertex(String kmer, boolean novelty) {
        this.kmer = kmer;
        this.isNovel = novelty;
    }

    public void setFlag(String flag) { flags.add(flag); }
    public void unsetFlag(String flag) { flags.remove(flag); }
    public boolean flagIsSet(String flag) { return flags.contains(flag); }

    public Set<String> getFlags() { return flags; }
    public void setFlags(Set<String> flags) { this.flags = flags; }

    public String getKmer() { return kmer; }

    public boolean isNovel() { return isNovel; }
    public void setNovel() { isNovel = true; }

    public boolean isAnnotated() { return isNovel || flags.size() > 0; }

    public void setMaternalLocations(Set<Interval> l) {
        ml.clear();
        ml.addAll(l);
    }

    public void setPaternalLocations(Set<Interval> l) {
        pl.clear();
        pl.addAll(l);
    }

    public Set<Interval> getMaternalLocations() { return ml; }
    public Set<Interval> getPaternalLocations() { return pl; }

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
