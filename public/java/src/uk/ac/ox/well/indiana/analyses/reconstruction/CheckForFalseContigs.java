package uk.ac.ox.well.indiana.analyses.reconstruction;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class CheckForFalseContigs extends Tool {
    @Argument(fullName="contigTable", shortName="ct", doc="Contig table")
    public File CONTIG_TABLE;

    @Argument(fullName="reference", shortName="R", doc="Reference sequence")
    public IndexedFastaSequenceFile REFERENCE;

    @Output
    public PrintStream out;

    private Set<String> loadReferenceSequence() {
        Set<String> ref = new HashSet<String>();
        ReferenceSequence seq;

        while ((seq = REFERENCE.nextSequence()) != null) {
            ref.add((new String(seq.getBases())).toUpperCase());
        }

        return ref;
    }

    @Override
    public void execute() {
        Set<String> ref = loadReferenceSequence();

        int found = 0, total = 0, id = 1;

        TableReader tr = new TableReader(CONTIG_TABLE);
        for (Map<String, String> te : tr) {
            String fw = te.get("contig");
            String rc = SequenceUtils.reverseComplement(fw);

            boolean isFound = false;

            for (String seq : ref) {
                if (seq.contains(fw) || seq.contains(rc)) {
                    found++;
                    isFound = true;
                    break;
                }
            }

            if (!isFound) {
                out.println(">missing_contig_" + id);
                out.println(fw);
            }

            total++;

            if (tr.size() > 10 && total % (tr.size() / 10) == 0) {
                log.info("Processed {}/{} records", total, tr.size());
            }

            id++;
        }

        log.info("Loaded {} contigs, found {} in reference", total, found);
    }
}
