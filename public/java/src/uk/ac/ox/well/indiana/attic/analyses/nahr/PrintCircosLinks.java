package uk.ac.ox.well.indiana.attic.analyses.nahr;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

import java.io.File;
import java.io.PrintStream;
import java.util.Map;

public class PrintCircosLinks extends Module {
    @Argument(fullName="ann", shortName="ann", doc="Contig annotations")
    public File ANN;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TableReader tr = new TableReader(ANN);
        for (Map<String, String> te : tr) {
            String contigName = te.get("contigName");

            String[] ref0ToCanonicalExact = (te.get("ref0ToCanonicalExact").equals("NA") || te.get("ref0ToCanonicalExact").equals("*:0-0") ? "NA:0-0" : te.get("ref0ToCanonicalExact")).split("[:-]");
            String[] ref1ToCanonicalExact = (te.get("ref1ToCanonicalExact").equals("NA") || te.get("ref1ToCanonicalExact").equals("*:0-0") ? "NA:0-0" : te.get("ref1ToCanonicalExact")).split("[:-]");

            if (ref0ToCanonicalExact[1].equals("") || ref0ToCanonicalExact[1].equals(" ")) { ref0ToCanonicalExact[1] = "0"; }
            if (ref0ToCanonicalExact[2].equals("") || ref0ToCanonicalExact[2].equals(" ")) { ref0ToCanonicalExact[2] = "0"; }
            if (ref1ToCanonicalExact[1].equals("") || ref1ToCanonicalExact[1].equals(" ")) { ref1ToCanonicalExact[1] = "0"; }
            if (ref1ToCanonicalExact[2].equals("") || ref1ToCanonicalExact[2].equals(" ")) { ref1ToCanonicalExact[2] = "0"; }

            out.println(te.get("sampleName") + "_" + te.get("accession") + "_" + contigName + " " + ref0ToCanonicalExact[0] + " " + ref0ToCanonicalExact[1] + " " + ref0ToCanonicalExact[2] + " radius1=0.8r,svgid=" + contigName);
            out.println(te.get("sampleName") + "_" + te.get("accession") + "_" + contigName + " " + ref1ToCanonicalExact[0] + " " + ref1ToCanonicalExact[1] + " " + ref1ToCanonicalExact[2] + " radius2=0.6r,svgid=" + contigName);
        }
    }
}
