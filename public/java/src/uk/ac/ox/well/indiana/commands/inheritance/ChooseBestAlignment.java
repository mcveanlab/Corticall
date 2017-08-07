package uk.ac.ox.well.indiana.commands.inheritance;

import com.google.common.base.Joiner;
import htsjdk.samtools.*;
import org.apache.commons.math3.util.MathUtils;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by kiran on 07/08/2017.
 */
public class ChooseBestAlignment extends Module {
    @Argument(fullName="sam", shortName="s", doc="SAM file")
    public ArrayList<SamReader> SAMS;

    @Output
    public PrintStream out;

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

                if (sr != null && sr.getCigar().getCigarElements().size() == 1) {
                    String[] pieces = sr.getReferenceName().replaceAll("_v3", "").split("_");
                    String chrName = "Pf3D7_" + pieces[pieces.length - 1] + "_v3";

                    //log.info("{} {}", contigName, sr.getSAMString());
                    out.println(Joiner.on(" ").join(chrName, sr.getStart(), sr.getEnd(), pieces[0]));
                }
            }
        }

        /*
chr - M76611 MT 0 5967 white
chr - PFC10_API_IRAB API 0 34242 white
chr - Pf3D7_01_v3 1 0 640851 white
chr - Pf3D7_02_v3 2 0 947102 white
chr - Pf3D7_03_v3 3 0 1067971 white
chr - Pf3D7_04_v3 4 0 1200490 white
chr - Pf3D7_05_v3 5 0 1343557 white
chr - Pf3D7_06_v3 6 0 1418242 white
chr - Pf3D7_07_v3 7 0 1445207 white
chr - Pf3D7_08_v3 8 0 1472805 white
chr - Pf3D7_09_v3 9 0 1541735 white
chr - Pf3D7_10_v3 10 0 1687656 white
chr - Pf3D7_11_v3 11 0 2038340 white
chr - Pf3D7_12_v3 12 0 2271494 white
chr - Pf3D7_13_v3 13 0 2925236 white
chr - Pf3D7_14_v3 14 0 3291936 white
chr - Pf3D7_00_v3 U 0 650000 white
         */
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
