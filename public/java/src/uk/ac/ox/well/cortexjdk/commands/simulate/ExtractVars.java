package uk.ac.ox.well.cortexjdk.commands.simulate;

import com.google.common.base.Joiner;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.gff.GFF3;
import uk.ac.ox.well.cortexjdk.utils.io.gff.GFF3Record;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by kiran on 26/10/2017.
 */
public class ExtractVars extends Module {
    @Argument(fullName="ref", shortName="r", doc="Reference")
    public KmerLookup REF;

    @Argument(fullName="gff", shortName="g", doc="GFF")
    public GFF3 GFF;

    @Argument(fullName="ids", shortName="i", doc="IDs")
    public HashSet<String> IDS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        List<Pair<String, String>> vars = new ArrayList<>();
        vars.addAll(loadVars(REF, GFF, IDS));

        for (Pair<String, String> p : vars) {
            out.println(">" + p.getFirst());
            out.println(p.getSecond());
        }
    }

    private List<Pair<String, String>> loadVars(KmerLookup kl, GFF3 gff, Set<String> ids) {
        Set<Pair<String, String>> vars = new HashSet<>();

        for (GFF3Record gr : gff) {
            if (ids.contains(gr.getAttribute("ID"))) {
                String seq = kl.getReferenceSequence().getSubsequenceAt(gr.getSeqid(), gr.getStart() - 500, gr.getEnd() + 500).getBaseString();

                if (gr.getStrand() == GFF3Record.Strand.NEGATIVE) {
                    seq = SequenceUtils.reverseComplement(seq);
                }

                vars.add(new Pair<>(gr.getAttribute("ID"), seq.toUpperCase()));
            }
        }

        return new ArrayList<>(vars);
    }
}
