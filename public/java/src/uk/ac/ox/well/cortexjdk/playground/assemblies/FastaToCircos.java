package uk.ac.ox.well.cortexjdk.playground.assemblies;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;

import java.io.PrintStream;

/**
 * Created by kiran on 01/07/2017.
 */
public class FastaToCircos extends Module {
    @Argument(fullName="sequence", shortName="s", doc="Sequences")
    public FastaSequenceFile SEQUENCE;

    @Argument(fullName="findPattern", shortName="fp", doc="Find regex")
    public String FIND_PATTERN;

    @Argument(fullName="replacePattern", shortName="rp", doc="Find regex")
    public String REPLACE_PATTERN;

    @Argument(fullName="breakAtNs", shortName="b", doc="Break at Ns")
    public Boolean BREAK_AT_NS = false;

    @Argument(fullName="keepNames", shortName="k", doc="Keep contig names")
    public Boolean KEEP_NAMES = false;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        ReferenceSequence rseq;
        while ((rseq = SEQUENCE.nextSequence()) != null) {
            String draftName = rseq.getName().split("\\s+")[0];
            String newName = KEEP_NAMES ? draftName : draftName.replaceAll(FIND_PATTERN, REPLACE_PATTERN);

            String[] fragments = new String[] { rseq.getBaseString() };
            if (BREAK_AT_NS) {
                fragments = rseq.getBaseString().split("[Nn]+");
            }

            int pos = 1;
            for (int i = 0; i < fragments.length; i++) {
                String fragment = fragments[i];
                out.println(Joiner.on(" ").join(newName, pos, pos + fragment.length(), draftName + (fragments.length == 1 ? "" : "_" + i)));

                pos += fragment.length() + 1;
            }
        }
    }
}
