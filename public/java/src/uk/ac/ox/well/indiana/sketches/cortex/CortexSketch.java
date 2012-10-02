package uk.ac.ox.well.indiana.sketches.cortex;

import uk.ac.ox.well.indiana.sketches.Sketch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;

import java.util.HashSet;

public class CortexSketch extends Sketch {
    @Argument(fullName="cortexGraph", shortName="cg", doc="A binary Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    private HashSet<CortexRecord> records = new HashSet<CortexRecord>();

    public void setup() {
        System.out.println(CORTEX_GRAPH);

        for (CortexRecord cr : CORTEX_GRAPH) {
            System.out.println(cr);
        }

        size(400, 400);
    }

    public void draw() {
        if (mousePressed) {
            ellipse(mouseX, mouseY, 20, 20);
        }
    }
}
