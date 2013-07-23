package uk.ac.ox.well.indiana.utils.assembly;

import com.google.common.base.Joiner;
import org.jgrapht.graph.DefaultEdge;

import java.util.Set;
import java.util.TreeSet;

public class CortexEdge extends DefaultEdge {
    private Set<Integer> colors = new TreeSet<Integer>();

    public void addColor(int color) {
        colors.add(color);
    }

    public String getLabel() {
        return Joiner.on(",").join(colors);
    }
}
