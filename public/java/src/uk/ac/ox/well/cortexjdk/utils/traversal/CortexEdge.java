package uk.ac.ox.well.cortexjdk.utils.traversal;

import org.jgrapht.graph.DefaultWeightedEdge;

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
}
