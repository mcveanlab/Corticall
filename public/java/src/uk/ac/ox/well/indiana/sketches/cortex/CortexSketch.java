package uk.ac.ox.well.indiana.sketches.cortex;

import uk.ac.ox.well.indiana.sketches.Sketch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;

public class CortexSketch extends Sketch {
    @Argument(fullName="cortexGraph", shortName="cg", doc="A binary Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    public void setup() {
        System.out.println(CORTEX_GRAPH);

        for (int i = 0; i < 10; i++) {
            CortexRecord cr = CORTEX_GRAPH.next();

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
