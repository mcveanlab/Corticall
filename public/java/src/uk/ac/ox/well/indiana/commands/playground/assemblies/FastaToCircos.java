package uk.ac.ox.well.indiana.commands.playground.assemblies;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

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
    public Boolean BREAK_AT_NS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        ReferenceSequence rseq;
        while ((rseq = SEQUENCE.nextSequence()) != null) {
            String draftName = rseq.getName().split("\\s+")[0];
            String newName = draftName.replaceAll(FIND_PATTERN, REPLACE_PATTERN);

            String[] fragments = new String[] { rseq.getBaseString() };
            if (BREAK_AT_NS) {
                fragments = rseq.getBaseString().split("[Nn]+");
            }

            int pos = 1;
            for (int i = 0; i < fragments.length; i++) {
                String fragment = fragments[i];
                out.println(Joiner.on(" ").join(newName, pos, pos + fragment.length(), draftName + "_" + i));

                pos += fragment.length();
            }
        }
    }
}
