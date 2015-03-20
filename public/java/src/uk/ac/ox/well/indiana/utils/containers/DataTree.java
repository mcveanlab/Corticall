package uk.ac.ox.well.indiana.utils.containers;

import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

public class DataTree {
    private Map<Object, DataTree> children;

    public void add(Object... objects) {
        add(objects, 0, objects.length);
    }

    public void add(Object[] objects, int start, int stop) {
        if (start < stop) {
            if (children == null) {
                children = new TreeMap<Object, DataTree>();
            }

            if (!children.containsKey(objects[start])) {
                children.put(objects[start], new DataTree());
            }

            children.get(objects[start]).add(objects, start + 1, stop);
        }
    }

    public DataTree get(Object o) {
        return children.get(o);
    }

    public Set<Object> get(Object... objects) {
        DataTree tree = this;

        for (int i = 0; i < objects.length; i++) {
            tree = tree.get(objects[i]);
        }

        return tree.children.keySet();
    }

    public boolean has(Object o) {
        return children.containsKey(o);
    }

    public boolean has(Object... objects) {
        DataTree tree = this;

        for (int i = 0; i < objects.length; i++) {
            if (!tree.has(objects[i])) {
                return false;
            }

            tree = tree.get(objects[i]);
        }

        return true;
    }
}
