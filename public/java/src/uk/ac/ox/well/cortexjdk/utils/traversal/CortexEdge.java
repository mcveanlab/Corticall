package uk.ac.ox.well.cortexjdk.utils.traversal;

import org.jgrapht.graph.DefaultWeightedEdge;

import java.util.*;

/**
 * Created by kiran on 10/05/2017.
 */
public class CortexEdge extends DefaultWeightedEdge {
    private Set<CortexVertex> vertices = new LinkedHashSet<>();

    private int color = -1;
    private double weight = 0.0;

    public CortexEdge() {}

    public CortexEdge(CortexVertex s, CortexVertex t, int color, double weight) {
        vertices.add(s);
        vertices.add(t);

        this.color = color;
        this.weight = weight;
    }

    public int getColor() { return color; }
    public void setColor(int color) { this.color = color; }

    @Override
    public double getWeight() { return weight; }
    public void setWeight(double weight) { this.weight = weight; }

    @Override
    public String toString() {
        return "CortexEdge{" +
                "color=" + color +
                ", weight=" + weight +
                ", vertices=" + vertices +
                '}';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        CortexEdge that = (CortexEdge) o;

        if (color != that.color) return false;
        if (Double.compare(that.weight, weight) != 0) return false;
        return vertices != null ? vertices.equals(that.vertices) : that.vertices == null;
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = vertices != null ? vertices.hashCode() : 0;
        result = 31 * result + color;
        temp = Double.doubleToLongBits(weight);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }
}
