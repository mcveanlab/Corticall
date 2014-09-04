package uk.ac.ox.well.indiana.attic.analyses.nahr;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
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
import java.util.*;

public class SelectContigsByMetrics extends Module {
    @Argument(fullName="contigMetrics", shortName="cm", doc="Contig metrics")
    public File ANN;

    @Argument(fullName="constraint", shortName="c", doc="A JEXL constraint to apply when selecting contigs to consider")
    public ArrayList<Expression> CONSTRAINTS;

    @Argument(fullName="alignedContigs", shortName="tc", doc="Aligned contigs (BAM)")
    public SAMFileReader CONTIGS;

    @Argument(fullName="sampleName", shortName="sn", doc="Sample name")
    public String SAMPLE_NAME;

    @Argument(fullName="accession", shortName="acc", doc="Accession")
    public String ACCESSION;

    @Argument(fullName="classifications", shortName="cl", doc="Classification files")
    public LinkedHashMap<String, File> CLASSIFICATIONS;

    @Output
    public PrintStream out;

    @Output(fullName="statsOut", shortName="so", doc="Stats out")
    public PrintStream sout;

    @Output(fullName="contigsOut", shortName="co", doc="Contig names out")
    public PrintStream cout;

    @Override
    public void execute() {
        /*
        Set<String> knownEvents = new HashSet<String>();
        for (SAMRecord read : TRAINING) {
            knownEvents.add(read.getReadName());
        }
        */

        Map<String, Interval> contigLoci = new HashMap<String, Interval>();
        for (SAMRecord contig : CONTIGS) {
            Interval it = new Interval(contig.getReferenceName(), contig.getAlignmentStart(), contig.getAlignmentEnd());

            contigLoci.put(contig.getReadName(), it);
        }

        Map<String, IntervalTreeMap<String>> its = new HashMap<String, IntervalTreeMap<String>>();
        its.put(ACCESSION, new IntervalTreeMap<String>());

        int knownEvents = 0;
        for (String type : CLASSIFICATIONS.keySet()) {
            TableReader ctr = new TableReader(CLASSIFICATIONS.get(type));

            for (Map<String, String> cte : ctr) {
                String[] samplePieces = cte.get("sample").split("/");
                String acc = samplePieces[2];
                String chrom = cte.get("chrom");
                Integer start = Integer.valueOf(cte.get("co_pos_min"));
                Integer stop = Integer.valueOf(cte.get("co_pos_max"));

                Interval it = new Interval(chrom, start, stop);

                if (!its.containsKey(acc)) {
                    its.put(acc, new IntervalTreeMap<String>());
                }

                if (acc.equals("all")) {
                    its.get(ACCESSION).put(it, type);
                }

                its.get(acc).put(it, type);

                if (acc.equals("all") || ACCESSION.equals(acc)) {
                    knownEvents++;
                }
            }
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
                /*
                te.put("isKnownAHR", "0");

                if (knownEvents.contains(te.get("contigName"))) {
                    knownContigs++;

                    te.put("isKnownAHR", "1");
                }
                */

                selectedContigs++;

                for (String type : CLASSIFICATIONS.keySet()) {
                    te.put(type, "0");
                }

                te.put("isUnclassified", "0");

                Interval it = contigLoci.get(te.get("contigName"));

                boolean isClassified = false;
                if (it != null && its.get(ACCESSION).containsOverlapping(it)) {
                    Set<String> types = new TreeSet<String>();
                    types.addAll(its.get(ACCESSION).getOverlapping(it));

                    for (String type : types) {
                        te.put(type, "1");
                        isClassified = true;
                    }

                    knownContigs++;
                }

                if (!isClassified) {
                    te.put("isUnclassified", "1");
                }

                tw.addEntry(te);

                cout.println(te.get("contigName"));
            }

            numContigs++;
        }

        log.info("Found {}/{} (~{}%) contigs that met criteria, {} training contigs recovered.", selectedContigs, numContigs, String.format("%.2f", 100.0f * (float) selectedContigs / (float) numContigs), knownContigs);

        sout.println(Joiner.on("\t").join(SAMPLE_NAME, numContigs, selectedContigs, knownEvents, knownContigs));
    }
}
