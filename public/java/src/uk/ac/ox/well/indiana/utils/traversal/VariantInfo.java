package uk.ac.ox.well.indiana.utils.traversal;

/**
 * Created by kiran on 31/08/2015.
 */
public class VariantInfo {
    public String variantId;

    public String vchr;
    public int vstart;
    public int vstop;

    public String type;
    public String denovo;
    public String nahr;
    public int gcindex;
    public String ref;
    public String alt;
    public String leftFlank;
    public String rightFlank;

    public boolean found = false;
    public boolean matches = false;

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        VariantInfo that = (VariantInfo) o;

        if (vstart != that.vstart) return false;
        if (vstop != that.vstop) return false;
        if (gcindex != that.gcindex) return false;
        if (variantId != null ? !variantId.equals(that.variantId) : that.variantId != null) return false;
        if (vchr != null ? !vchr.equals(that.vchr) : that.vchr != null) return false;
        if (type != null ? !type.equals(that.type) : that.type != null) return false;
        if (denovo != null ? !denovo.equals(that.denovo) : that.denovo != null) return false;
        if (nahr != null ? !nahr.equals(that.nahr) : that.nahr != null) return false;
        if (ref != null ? !ref.equals(that.ref) : that.ref != null) return false;
        if (alt != null ? !alt.equals(that.alt) : that.alt != null) return false;
        if (leftFlank != null ? !leftFlank.equals(that.leftFlank) : that.leftFlank != null) return false;
        return !(rightFlank != null ? !rightFlank.equals(that.rightFlank) : that.rightFlank != null);
    }

    @Override
    public int hashCode() {
        int result = variantId != null ? variantId.hashCode() : 0;
        result = 31 * result + (vchr != null ? vchr.hashCode() : 0);
        result = 31 * result + vstart;
        result = 31 * result + vstop;
        result = 31 * result + (type != null ? type.hashCode() : 0);
        result = 31 * result + (denovo != null ? denovo.hashCode() : 0);
        result = 31 * result + (nahr != null ? nahr.hashCode() : 0);
        result = 31 * result + gcindex;
        result = 31 * result + (ref != null ? ref.hashCode() : 0);
        result = 31 * result + (alt != null ? alt.hashCode() : 0);
        result = 31 * result + (leftFlank != null ? leftFlank.hashCode() : 0);
        result = 31 * result + (rightFlank != null ? rightFlank.hashCode() : 0);
        return result;
    }

    @Override
    public String toString() {
        return "VariantInfo{" +
                "variantId='" + variantId + '\'' +
                ", vchr='" + vchr + '\'' +
                ", vstart=" + vstart +
                ", vstop=" + vstop +
                ", type='" + type + '\'' +
                ", denovo='" + denovo + '\'' +
                ", nahr='" + nahr + '\'' +
                ", gcindex=" + gcindex +
                ", ref='" + ref + '\'' +
                ", alt='" + alt + '\'' +
                ", leftFlank='" + leftFlank + '\'' +
                ", rightFlank='" + rightFlank + '\'' +
                '}';
    }
}

