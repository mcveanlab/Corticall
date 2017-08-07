package uk.ac.ox.well.indiana.commands.inheritance;

import htsjdk.samtools.*;
import org.apache.commons.math3.util.MathUtils;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.File;
import java.util.*;

/**
 * Created by kiran on 07/08/2017.
 */
public class ChooseBestAlignment extends Module {
    @Argument(fullName="sam", shortName="s", doc="SAM file")
    public ArrayList<SamReader> SAMS;

    //@Output
    //public File out;

    @Override
    public void execute() {
        List<Map<String, List<SAMRecord>>> contigs = new ArrayList<>();
        Set<String> chosenContigs = new TreeSet<>();

        for (int i = 0; i < SAMS.size(); i++) {
            contigs.add(new HashMap<>());

            for (SAMRecord sr : SAMS.get(i)) {
                if (!contigs.get(i).containsKey(sr.getReadName())) {
                    contigs.get(i).put(sr.getReadName(), new ArrayList<>());
                }

                contigs.get(i).get(sr.getReadName()).add(sr);
                chosenContigs.add(sr.getReadName());
            }
        }

        for (String contigName : chosenContigs) {
            if (contigs.get(0).containsKey(contigName) && contigs.get(1).containsKey(contigName) &&
                contigs.get(0).get(contigName).size() == 1 && contigs.get(1).get(contigName).size() == 1) {
                SAMRecord sr0 = contigs.get(0).get(contigName).get(0);
                SAMRecord sr1 = contigs.get(1).get(contigName).get(0);

                SAMRecord sr = chooseBetterAlignment(sr0, sr1);

                if (sr != null) {
                    log.info("{} {}", contigName, sr.getSAMString());
                }
            }
        }
    }

    private SAMRecord chooseBetterAlignment(SAMRecord s0, SAMRecord s1) {
        int l0 = 0, l1 = 0;

        for (CigarElement ce : s0.getCigar().getCigarElements()) {
            if (ce.getOperator().equals(CigarOperator.MATCH_OR_MISMATCH)) {
                l0 += ce.getLength();
            } else {
                l0 -= ce.getLength();
            }
        }

        for (CigarElement ce : s1.getCigar().getCigarElements()) {
            if (ce.getOperator().equals(CigarOperator.MATCH_OR_MISMATCH)) {
                l1 += ce.getLength();
            } else {
                l1 -= ce.getLength();
            }
        }

        if (s0.getCigar().equals(s1.getCigar())) {
            return null;
        }

        return l0 > l1 ? s0 : s1;
    }
}
