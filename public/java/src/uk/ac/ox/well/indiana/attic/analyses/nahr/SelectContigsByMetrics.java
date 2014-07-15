package uk.ac.ox.well.indiana.attic.analyses.nahr;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
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
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class SelectContigsByMetrics extends Module {
    @Argument(fullName="contigMetrics", shortName="cm", doc="Contig metrics")
    public File ANN;

    @Argument(fullName="constraint", shortName="c", doc="A JEXL constraint to apply when selecting contigs to consider")
    public ArrayList<Expression> CONSTRAINTS;

    @Argument(fullName="trainingContigs", shortName="tc", doc="Training contigs (BAM)")
    public SAMFileReader TRAINING;

    @Argument(fullName="sampleName", shortName="sn", doc="Sample name")
    public String SAMPLE_NAME;

    @Output
    public PrintStream out;

    @Output(fullName="statsOut", shortName="so", doc="Stats out")
    public PrintStream sout;

    @Output(fullName="contigsOut", shortName="co", doc="Contig names out")
    public PrintStream cout;

    @Override
    public void execute() {
        Set<String> knownEvents = new HashSet<String>();
        for (SAMRecord read : TRAINING) {
            knownEvents.add(read.getReadName());
        }

        TableReader tr = new TableReader(ANN);
        TableWriter tw = new TableWriter(out);

        int numContigs = 0;
        int selectedContigs = 0;
        int knownContigs = 0;

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
                te.put("isKnownAHR", "0");

                if (knownEvents.contains(te.get("contigName"))) {
                    knownContigs++;

                    te.put("isKnownAHR", "1");
                }
                selectedContigs++;

                tw.addEntry(te);

                cout.println(te.get("contigName"));
            }

            numContigs++;
        }

        log.info("Found {}/{} (~{}%) contigs that met criteria, {} training contigs recovered.", selectedContigs, numContigs, String.format("%.2f", 100.0f * (float) selectedContigs / (float) numContigs), knownContigs);

        sout.println(Joiner.on("\t").join(SAMPLE_NAME, numContigs, selectedContigs, knownEvents.size(), knownContigs));
    }
}
