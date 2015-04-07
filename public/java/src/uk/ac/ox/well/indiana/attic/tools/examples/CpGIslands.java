package uk.ac.ox.well.indiana.attic.tools.examples;

import com.google.common.base.Joiner;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.math.MoreMathUtils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CpGIslands extends Module {
    @Argument(fullName="sequence", shortName="s", doc="Sequence to classify")
    public String SEQUENCE;

    @Output
    public PrintStream out;

    private void printMatrix(double[][] matrix, String[] rowNames, String[] colNames) {
        out.println("\t" + Joiner.on("\t").join(colNames));

        for (int r = 0; r < rowNames.length; r++) {
            List<Double> values = new ArrayList<Double>();

            for (int c = 0; c < colNames.length; c++) {
                values.add(matrix[r][c]);
            }

            out.println(rowNames[r] + "\t" + Joiner.on("\t").join(values));
        }
    }

    private void viterbi(String sequence, double[][] transitionmatrix, double[][] emissionmatrix, String[] states) {
        double[][] v = makeViterbimat(sequence, transitionmatrix, emissionmatrix);

        int[] mostprobablestatepath = new int[v.length];
        for (int r = 0; r < v.length; r++) {
            mostprobablestatepath[r] = (v[r][0] > v[r][1]) ? 0 : 1;
        }

        char prevnucleotide = sequence.charAt(0);
        int prevmostprobablestate = mostprobablestatepath[0];
        String prevmostprobablestatename = states[prevmostprobablestate];
        int startpos = 0;

        for (int i = 1; i < sequence.length(); i++) {
            char nucleotide = sequence.charAt(i);
            int mostprobablestate = mostprobablestatepath[i];
            String mostprobablestatename = states[mostprobablestate];
            if (!mostprobablestatename.equals(prevmostprobablestatename)) {
                log.info("Positions {}-{} Most probable state = {}", startpos, i - 1, prevmostprobablestatename);
                startpos = i;
            }
            prevnucleotide = nucleotide;
            prevmostprobablestatename = mostprobablestatename;
        }

        log.info("Positions {}-{} Most probable state = {}", startpos, sequence.length() - 1, prevmostprobablestatename);
    }

    private double[] innerProduct(double[] a, double[] b) {
        double[] c = new double[a.length];

        for (int i = 0; i < a.length; i++) {
            c[i] = a[i] * b[i];
        }

        return c;
    }

    private double[][] makeViterbimat(String sequence, double[][] transitionmatrix, double[][] emissionmatrix) {
        int numstates = transitionmatrix.length;

        double[][] v = new double[sequence.length()][numstates];
        v[0][0] = 1;

        for (int i = 1; i < sequence.length(); i++) {
            for (int l = 0; l < numstates; l++) {
                int emissionIndex = 0;
                switch (sequence.charAt(i)) {
                    case 'A': emissionIndex = 0; break;
                    case 'C': emissionIndex = 1; break;
                    case 'G': emissionIndex = 2; break;
                    case 'T': emissionIndex = 3; break;
                    default: emissionIndex = 0; break;
                }

                double[] transitionmatrixcol = new double[transitionmatrix[0].length];
                for (int r = 0; r < transitionmatrix.length; r++) {
                    transitionmatrixcol[r] = transitionmatrix[r][l];
                }

                double statelprobnucleotidei = emissionmatrix[l][emissionIndex];
                v[i][l] = statelprobnucleotidei * MoreMathUtils.max(innerProduct(v[i-1], transitionmatrixcol));
            }
        }

        return v;
    }

    @Override
    public void execute() {
        // Define states
        String[] states = new String[] { "AT-rich", "GC-rich" };

        // Define emissions
        String[] nucleotides = new String[] { "A", "C", "G", "T" };

        // Define transition matrix
        double[] ATrichprobs = new double[] { 0.7, 0.3 };
        double[] GCrichprobs = new double[] { 0.1, 0.9 };

        double[][] thetransitionmatrix = new double[2][2];
        thetransitionmatrix[0] = ATrichprobs;
        thetransitionmatrix[1] = GCrichprobs;

        // Define emission matrix
        double[] ATrichstateprobs = new double[] { 0.39, 0.1, 0.1, 0.41 };
        double[] GCrichstateprobs = new double[] { 0.1, 0.41, 0.39, 0.1 };

        double[][] theemissionmatrix = new double[2][4];
        theemissionmatrix[0] = ATrichstateprobs;
        theemissionmatrix[1] = GCrichstateprobs;

        //printMatrix(thetransitionmatrix, states, states);
        //printMatrix(theemissionmatrix, states, nucleotides);

        //double[][] v = makeViterbimat(SEQUENCE, thetransitionmatrix, theemissionmatrix);
        viterbi(SEQUENCE, thetransitionmatrix, theemissionmatrix, states);

        //log.info("{}: {}", SEQUENCE, SEQUENCE.split("").length);
        //log.info("{} {}", v.length, v[0].length);
        //log.info("{}", Joiner.on(" ").join(SEQUENCE.split("")));

        //printMatrix(v, Arrays.copyOfRange(SEQUENCE.split(""), 1, SEQUENCE.length() + 1), states);
    }
}
