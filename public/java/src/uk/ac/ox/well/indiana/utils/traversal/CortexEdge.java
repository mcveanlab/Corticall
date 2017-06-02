package uk.ac.ox.well.indiana.utils.traversal;

import org.jgrapht.graph.DefaultWeightedEdge;

import java.util.BitSet;

/**
 * Created by kiran on 10/05/2017.
 */
public class CortexEdge extends DefaultWeightedEdge {
    private int color = -1;
    private double weight = 0.0;

    public CortexEdge() {}

    public CortexEdge(int color, double weight) {
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
                '}';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        CortexEdge that = (CortexEdge) o;

        if (color != that.color) return false;
        return Double.compare(that.weight, weight) == 0;

    }

    /*
    @Override
    public int hashCode() {
        int result;
        long temp;
        result = color;
        temp = Double.doubleToLongBits(weight);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }
    */
}
