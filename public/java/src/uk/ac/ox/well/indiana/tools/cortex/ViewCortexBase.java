package uk.ac.ox.well.indiana.tools.cortex;

import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlContext;
import org.apache.commons.jexl2.MapContext;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;

import java.io.PrintStream;
import java.util.ArrayList;

public abstract class ViewCortexBase extends Tool {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="constraint", shortName="c", doc="A JEXL constraint to apply when selecting kmers to consider")
    public ArrayList<Expression> CONSTRAINTS;

    @Output
    public PrintStream out;

    public boolean satisfiesConstraints(CortexRecord cr) {
        JexlContext jexlContext = new MapContext();

        for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
            String sampleName = CORTEX_GRAPH.getColor(color).getSampleName();
            int coverage = cr.getCoverages()[color];

            jexlContext.set("color." + color, coverage);
            jexlContext.set(sampleName, coverage);
        }

        boolean allConstraintsSatisfied = true;
        for (Expression e : CONSTRAINTS) {
            Boolean constraintSatisfied = (Boolean) e.evaluate(jexlContext);

            if (!constraintSatisfied) {
                allConstraintsSatisfied = false;
                break;
            }
        }

        return allConstraintsSatisfied;
    }
}
