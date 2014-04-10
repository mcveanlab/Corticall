package uk.ac.ox.well.indiana.attic.alignment.pairwise;

import com.google.common.base.Joiner;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static uk.ac.ox.well.indiana.utils.math.MoreMathUtils.max;

public class PairwiseAligner extends Module {
    @Argument(fullName="query", shortName="q", doc="Query sequence(s)")
    public ArrayList<String> QUERIES;

    @Argument(fullName="target", shortName="t", doc="Target sequence(s)")
    public ArrayList<String> TARGETS;

    @Argument(fullName="training", shortName="train", doc="Training data in BAM format", required=false)
    public SAMFileReader TRAINING;

    @Argument(fullName="match", shortName="epm", doc="Emission probability for a match")
    public Double MATCH = 0.9;

    @Argument(fullName="insertion", shortName="epi", doc="Emission probability for an insertion")
    public Double INSERTION = 0.05;

    @Argument(fullName="deletion", shortName="epd", doc="Emission probability for a deletion")
    public Double DELETION = 0.05;

    @Argument(fullName="delta", shortName="tpd", doc="Transition probability for first gap")
    public Double DELTA = 0.2;

    @Argument(fullName="epsilon", shortName="tpe", doc="Transition probability for remaining in gap")
    public Double EPSILON = 0.2;

    @Argument(fullName="tau", shortName="tpt", doc="Transition probability for end state")
    public Double TAU = 0.1;

    @Output
    public PrintStream out;

    private void printMatrix(double a[][], String target, String query, String label) {
        out.println(label + ":");

        for (int j = 0; j < target.length(); j++) { out.print("\t" + target.charAt(j)); } out.println();

        for (int i = 0; i < query.length(); i++) {
            out.print(query.charAt(i));
            for (int j = 0; j < target.length(); j++) {
                out.print("\t" + a[i][j]);
            }
            out.println();
        }
    }

    @Override
    public void execute() {
        if (TRAINING != null) {
            int matchesSeen = 0;
            int mismatchesSeen = 0;
            int deletionsSeen = 0;
            int insertionsSeen = 0;
            int matchLengthSum = 0;
            int deletionLengthSum = 0;
            int insertionLengthSum = 0;

            Pattern p = Pattern.compile("[0-9]+|[A-Z]|\\^[A-Z]+");

            for (SAMRecord read : TRAINING) {
                String md = read.getStringAttribute("MD");

                int numMatches = 0;
                int numMismatches = 0;
                int numDeletions = 0;
                int numInsertions = 0;

                for (CigarElement ce : read.getCigar().getCigarElements()) {
                    if (ce.getOperator().equals(CigarOperator.I)) {
                        numInsertions++;
                    } else if (ce.getOperator().equals(CigarOperator.D)) {
                        numDeletions++;
                    }
                }

                List<String> pieces = new ArrayList<String>();

                if (md != null) {
                    Matcher m = p.matcher(md);

                    while (m.find()) {
                        String piece = md.substring(m.start(), m.end());

                        if (piece.matches("[A-Z]") && !piece.contains("^")) {
                            numMismatches++;
                        } else {
                            numMatches += Integer.parseInt(piece);
                        }

                        pieces.add(piece);
                    }
                }

                String joined = Joiner.on(",").join(pieces);

                out.println(Joiner.on("\t").useForNull("").join(read.getCigarString(), md, joined, numMatches, numMismatches, numDeletions, numInsertions));
            }
        }

        for (String query : QUERIES) {
            for (String target : TARGETS) {
                double vM[][] = new double[query.length()][target.length()];
                double vI[][] = new double[query.length()][target.length()];
                double vD[][] = new double[query.length()][target.length()];

                vM[0][0] = 1.0;

                for (int i = 0; i < query.length(); i++) {
                    for (int j = 0; j < target.length(); j++) {
                        if (!(i == 0 && j == 0)) {
                            double vM_M = 0.0;
                            double vM_I = 0.0;
                            double vM_D = 0.0;

                            double vI_M = 0.0;
                            double vI_I = 0.0;

                            double vD_M = 0.0;
                            double vD_D = 0.0;

                            if (i - 1 >= 0 && j - 1 >= 0) {
                                vM_M = (1 - 2*DELTA - TAU)*vM[i - 1][j - 1];
                                vM_I = (1 - EPSILON - TAU)*vI[i - 1][j - 1];
                                vM_D = (1 - EPSILON - TAU)*vD[i - 1][j - 1];
                            }

                            if (i - 1 >= 0) {
                                vI_M = DELTA*vM[i - 1][j];
                                vI_I = EPSILON*vI[i - 1][j];
                            }

                            if (j - 1 >= 0) {
                                vD_M = DELTA*vM[i][j - 1];
                                vD_D = EPSILON*vD[i][j - 1];
                            }

                            vM[i][j] = MATCH*max(vM_M, vM_I, vM_D);
                            vI[i][j] = INSERTION*max(vI_M, vI_I);
                            vI[i][j] = DELETION*max(vD_M, vD_D);
                        }
                    }
                }

                printMatrix(vM, target, query, "M");
                printMatrix(vI, target, query, "I");
                printMatrix(vD, target, query, "D");
            }
        }
    }
}
