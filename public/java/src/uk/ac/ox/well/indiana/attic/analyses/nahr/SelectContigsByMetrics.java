package uk.ac.ox.well.indiana.attic.analyses.nahr;

import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlContext;
import org.apache.commons.jexl2.MapContext;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Map;

public class SelectContigsByMetrics extends Module {
    @Argument(fullName="contigMetrics", shortName="cm", doc="Contig metrics")
    public File ANN;

    @Argument(fullName="constraint", shortName="c", doc="A JEXL constraint to apply when selecting contigs to consider")
    public ArrayList<Expression> CONSTRAINTS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TableReader tr = new TableReader(ANN);
        TableWriter tw = new TableWriter(out);

        int numContigs = 0;
        int selectedContigs = 0;

        for (Map<String, String> te : tr) {
            JexlContext jexlContext = new MapContext();

            for (String key : te.keySet()) {
                jexlContext.set(key, te.get(key));
            }

            boolean allSatisfied = true;
            for (Expression e : CONSTRAINTS) {
                try {
                    Boolean satisfied = (Boolean) e.evaluate(jexlContext);

                    if (!satisfied) {
                        allSatisfied = false;
                        break;
                    }
                } catch (ClassCastException ex) {
                    log.error("{}", te);
                    throw new IndianaException("Problem evaluating JEXL expression for expression " + e + ": ", ex);
                }
            }

            if (allSatisfied) {
                tw.addEntry(te);

                selectedContigs++;
            }

            numContigs++;
        }

        log.info("Found {}/{} (~{}%) contigs that met criteria.", selectedContigs, numContigs, String.format("%.2f", 100.0f * (float) selectedContigs / (float) numContigs));
    }
}
