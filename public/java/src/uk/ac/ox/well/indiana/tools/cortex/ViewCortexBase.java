package uk.ac.ox.well.indiana.tools.cortex;

import net.sf.picard.reference.FastaSequenceFile;
import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlContext;
import org.apache.commons.jexl2.MapContext;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Map;

public abstract class ViewCortexBase extends Tool {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="constraint", shortName="c", doc="A JEXL constraint to apply when selecting kmers to consider")
    public ArrayList<Expression> CONSTRAINTS;

    @Argument(fullName="targetsOfInterest", shortName="toi", doc="One or more fasta file of targets to find", required=false)
    public ArrayList<FastaSequenceFile> TARGETS_OF_INTEREST;

    @Output
    public PrintStream out;

    protected Map<Integer, String> kmerMap = null;

    public boolean satisfiesConstraints(CortexRecord cr) {
        if (TARGETS_OF_INTEREST != null && kmerMap == null) {
            kmerMap = SequenceUtils.loadSequenceCodesAsAlphanumericallyLowestKmers(TARGETS_OF_INTEREST, CORTEX_GRAPH.getKmerSize());
        }

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

        if (allConstraintsSatisfied && kmerMap != null && !kmerMap.containsKey(cr.getKmerString().hashCode())) {
            allConstraintsSatisfied = false;
        }

        return allConstraintsSatisfied;
    }

    public String getKmerHomeContigName(CortexRecord cr) {
        if (kmerMap == null) {
            return null;
        }

        String homeName = kmerMap.get(cr.getKmerString().hashCode());

        if (homeName == null) {
            return "none";
        }

        return homeName;
    }
}
