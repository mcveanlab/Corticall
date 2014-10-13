package uk.ac.ox.well.indiana.commands.cortex;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.io.PrintStream;
import java.util.ArrayList;

public class print extends Module {
    @Argument(fullName="graph", shortName="g", doc="Cortex graph")
    public CortexGraph GRAPH;

    @Argument(fullName="record", shortName="r", doc="Record to find and print", required=false)
    public ArrayList<String> RECORDS;

    @Argument(fullName="headerOnly", shortName="H", doc="Only print the file header", required=false)
    public Boolean HEADER_ONLY = false;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        if (HEADER_ONLY) {
            out.println(GRAPH.getHeader());
        } else {
            if (RECORDS != null && !RECORDS.isEmpty()) {
                for (String kmer : RECORDS) {
                    CortexKmer ck = new CortexKmer(kmer);
                    CortexRecord cr = GRAPH.findRecord(ck);

                    if (cr != null) {
                        out.println(kmer + ": " + cr);
                    } else {
                        out.println(kmer + ": missing");
                    }
                }
            } else {
                for (CortexRecord cr : GRAPH) {
                    out.println(cr);
                }
            }
        }
    }
}