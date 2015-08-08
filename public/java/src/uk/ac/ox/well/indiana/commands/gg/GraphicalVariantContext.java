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

    public GraphicalVariantContext attribute(int c, String k, Object v) {
        attrs.get(c).put(k, v);

        return this;
    }

    public boolean hasAttribute(int c, String k) {
        return attrs.get(c).containsKey(k);
    }

    @SuppressWarnings("unchecked")
    public Object getAttribute(int c, String k) {
        /*
        if (c >= attrs.size()) {
            throw new IndianaException("Color '" + c + "' not found");
        }

        if (!attrs.get(c).containsKey(k)) {
            throw new IndianaException("Key '" + k + "' in c '" + c + "' not found");
        }
        */

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
}
