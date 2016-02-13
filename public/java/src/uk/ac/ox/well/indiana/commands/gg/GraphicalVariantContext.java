package uk.ac.ox.well.indiana.commands.gg;

import htsjdk.variant.variantcontext.VariantContext;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class GraphicalVariantContext {
    private List<Map<String, Object>> attrs = new ArrayList<Map<String, Object>>();

    public GraphicalVariantContext() {
        for (int i = 0; i <= 2; i++) {
            attrs.add(new TreeMap<String, Object>());
        }
    }

    public GraphicalVariantContext(GraphicalVariantContext g) {
        for (int i = 0; i <= 2; i++) {
            attrs.add(new TreeMap<String, Object>());

            addAttributes(i, g.getAttributes(i));
        }
    }

    public GraphicalVariantContext attribute(int c, String k, Object v) {
        attrs.get(c).put(k, v);

        return this;
    }

    public boolean hasAttribute(int c, String k) {
        return attrs.get(c).containsKey(k);
    }

    @SuppressWarnings("unchecked")
    public Object getAttribute(int c, String k) {
        return hasAttribute(c, k) ? attrs.get(c).get(k) : null;
    }

    public int getAttributeAsInt(int c, String k) {
        return hasAttribute(c, k) ? (Integer) getAttribute(c, k) : -1;
    }

    public String getAttributeAsString(int c, String k) {
        return hasAttribute(c, k) ? (String) getAttribute(c, k) : "";
    }

    public Boolean getAttributeAsBoolean(int c, String k) {
        return hasAttribute(c, k) ? (Boolean) getAttribute(c, k) : false;
    }

    public GraphicalVariantContext add(GraphicalVariantContext gvc) {
        for (int c = 0; c < gvc.attrs.size(); c++) {
            attrs.get(c).putAll(gvc.attrs.get(c));
        }

        return this;
    }

    public Map<String, Object> getAttributes(int c) {
        return attrs.get(c);
    }

    public void addAttributes(int c, Map<String, Object> attr) {
        attrs.get(c).putAll(attr);
    }

    public VariantContext.Type getType(int c) {
        String parentalAllele = getAttributeAsString(c, "parentalAllele");
        String childAllele = getAttributeAsString(c, "childAllele");

        if (parentalAllele.length() == 1 && childAllele.length() == 1) {
            return VariantContext.Type.SNP;
        } else if (parentalAllele.length() == 1 && childAllele.length() > 1) {
            return VariantContext.Type.INDEL;
        } else if (parentalAllele.length() > 1 && childAllele.length() == 1) {
            return VariantContext.Type.INDEL;
        } else if (parentalAllele.length() > 1 && childAllele.length() > 1) {
            return VariantContext.Type.MNP;
        }

        return VariantContext.Type.NO_VARIATION;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();

        sb.append(this.hashCode()).append("\n");

        for (int c = 0; c < 3; c++) {
            Map<String, Object> attr = attrs.get(c);

            for (String k : attr.keySet()) {
                sb.append(c).append(": ").append(k).append(" => ").append(SequenceUtils.truncate(attr.get(k).toString(), 100)).append("\n");
            }
        }

        return sb.toString();
    }

    /*
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        GraphicalVariantContext that = (GraphicalVariantContext) o;

        return !(attrs != null ? !attrs.equals(that.attrs) : that.attrs != null);

    }

    @Override
    public int hashCode() {
//        int result = kmer != null ? kmer.hashCode() : 0;
//        result = 31 * result + (isNovel ? 1 : 0);
//        return result;

        int result =               (attrs.get(0).containsKey("childAllele")         ? attrs.get(0).get("childAllele").hashCode()         : 0);
        result     = 31 * result + (attrs.get(0).containsKey("parentalAllele")      ? attrs.get(0).get("parentalAllele").hashCode()      : 0);
        result     = 31 * result + (attrs.get(0).containsKey("event")               ? attrs.get(0).get("event").hashCode()               : 0);
        result     = 31 * result + (attrs.get(0).containsKey("haplotypeBackground") ? attrs.get(0).get("haplotypeBackground").hashCode() : 0);
        result     = 31 * result + (attrs.get(0).containsKey("traversalStatus")     ? attrs.get(0).get("traversalStatus").hashCode()     : 0);


//        0 = {TreeMap$Entry@1548} "childAllele" -> "AAT"
//        1 = {TreeMap$Entry@1549} "childStretch" -> "CATTTATTCCTATCCACATTTCAAAGGTAAAACCTTCATTGTTATCAAATTTGTTAGGTTTTTGAAAATTTCTTTTTTTATAAGTATACTTTTTT"
//        2 = {TreeMap$Entry@1550} "event" -> "INS"
//        3 = {TreeMap$Entry@1551} "haplotypeBackground" -> "0"
//        4 = {TreeMap$Entry@1552} "novelKmersTotal" -> "1430"
//        5 = {TreeMap$Entry@1553} "novelKmersUsed" -> "47"
//        6 = {TreeMap$Entry@1554} "novelStretchAlignment" -> " size = 6"
//        7 = {TreeMap$Entry@1555} "parentalAllele" -> "T"
//        8 = {TreeMap$Entry@1556} "parentalPathAlignment" -> " size = 1"
//        9 = {TreeMap$Entry@1557} "parentalStretch" -> "CATTTATTCCTATCCACATTTCAAAGGTAAAACCTTCATTGTTATCATTTGTTAGGTTTTTGAAAATTTCTTTTTTTATAAGTATACTTTTTT"
//        10 = {TreeMap$Entry@1558} "score" -> "2"
//        11 = {TreeMap$Entry@1559} "start" -> "47"
//        12 = {TreeMap$Entry@1560} "stop" -> "47"
//        13 = {TreeMap$Entry@1561} "stretch" -> "ATCACCACAATCAATTTGATTATCACCACAATCAATTTGATTATCATCACAATCAATTTTATCATCACCATTTTTATTATTTTCCTTATTTTCTATCAATAATTCGTCACCCACAAACGGTACCGATGTTGACACTCCCTCAGGATATGCACCCTTCAATATATCATTTGTCTGTCCAAATATATTGACCTCGTCACTATTTACATCATACTGAATAACCATCGAAGGGTCATCAAACAATGATTCATATTTTTCCTCTGATATTAAGGAAAAGTCTACAGAGGAACTTTTTAAAAAACCATTAAATTTGACGGTTGACGAAAAAATATTAAAATTACTTGTTGTATATATATAATATCCCCACTCAGGAGAAAAATATGATATGCACTTTTCAAAATTGAATTGATAATTATTTATTTGAGCAAATTTATGACTCTTATAACATGTTGAATATATAAGAAAGTTGTTTTTTGAAAGCTCATATATCTGGTTCAATAACATTTTATTATTATGTTGCACATTTTCTTTTTTTTTTCTACCTTTTGTTGTGAAGGATTCATATTTGCCCTTATTTTCATAAGGAAATAACAAGAAATAAAATTCACACATCTTTTGAACATAAAATAAACATCCTCTCAAAAACAAACTTTTTATAAATTTCTTTAGAACATATTGTAATAATTTCGTATTTACTACAATTCGACATAATATAGAATATAATATTAATATAAGTGATGAATTCTCAATATATCTATAATGTTTATATTCATTAATCATTAATAAGATTTTACATATGGTTATATTCCTTACATTCGTTTTTAATACATTCTTTAATTCATCCATAGTTATTTCGGAACTATCTAAATTTAAATAGAAATCAATATATTTGTGATTATCTTTTAATTTCTTATATTTTAAAATACAATAACTTCTAACAAAATTTTTTATATAATCACTATTACGTATACATACTTTTAAAGCTTTTG
//        14 = {TreeMap$Entry@1562} "stretchLength" -> "4724"
//        15 = {TreeMap$Entry@1563} "stretchNum" -> "1"
//        16 = {TreeMap$Entry@1564} "traversalStatus" -> "complete"
    }
    */
}
