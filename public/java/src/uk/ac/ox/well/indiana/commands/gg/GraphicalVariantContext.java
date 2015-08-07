package uk.ac.ox.well.indiana.commands.gg;

import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;

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

    public Object getAttribute(int c, String k) {
        if (c >= attrs.size()) {
            throw new IndianaException("Color '" + c + "' not found");
        }

        if (!attrs.get(c).containsKey(k)) {
            throw new IndianaException("Key '" + k + "' in c '" + c + "' not found");
        }

        return attrs.get(c).get(k);
    }

    public int getAttributeAsInt(int c, String k) {
        return (Integer) getAttribute(c, k);
    }

    public String getAttributeAsString(int c, String k) {
        return (String) getAttribute(c, k);
    }

    public GraphicalVariantContext add(GraphicalVariantContext gvc) {
        for (int c = 0; c < gvc.attrs.size(); c++) {
            attrs.get(c).putAll(gvc.attrs.get(c));
        }

        return this;
    }
}
