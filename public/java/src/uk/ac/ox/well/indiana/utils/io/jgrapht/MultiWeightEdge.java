package uk.ac.ox.well.indiana.utils.io.jgrapht;

import org.jgrapht.graph.DefaultEdge;

import java.util.HashMap;
import java.util.Map;

public class MultiWeightEdge extends DefaultEdge {
    private Map<String, Integer> weights = new HashMap<String, Integer>();

    public void incrementWeight(String color) {
        if (weights.containsKey(color)) {
            weights.put(color, weights.get(color) + 1);
        } else {
            weights.put(color, 1);
        }
    }

    public Map<String, Integer> getWeights() {
        return weights;
    }

    public int totalWeight() {
        int weight = 0;

        for (String color : weights.keySet()) {
            weight += weights.get(color);
        }

        return weight;
    }
}
