package uk.ac.ox.well.indiana.utils.alignment.pairwise;

import com.google.common.base.Joiner;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.indiana.utils.math.MoreMathUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GlobalAligner {
    private double[][] tm;
    private double[][] em_match;
    private double em_indel;

    //private double[][] em;

    public GlobalAligner() { initialize(0.2, 0.1, 0.1); }

    public GlobalAligner(double delta, double epsilon, double tau) { initialize(delta, epsilon, tau); }

    private void initialize(double delta, double epsilon, double tau) {
        // Define transition matrix
        tm = new double[5][5];
        tm[0] = new double[] { 0, 1 - (2*delta) - tau, delta,   delta,   tau };
        tm[1] = new double[] { 0, 1 - (2*delta) - tau, delta,   delta,   tau };
        tm[2] = new double[] { 0, 1 - epsilon - tau,   epsilon, 0,       tau };
        tm[3] = new double[] { 0, 1 - epsilon - tau,   0,       epsilon, tau };
        tm[4] = new double[] { 0, 0,                   0,       0,       0   };

        // Define emission for matches
        em_match = new double[4][4];
        em_match[0] = new double[] { 0.50, 0.05, 0.15, 0.30 };
        em_match[1] = new double[] { 0.05, 0.50, 0.30, 0.15 };
        em_match[2] = new double[] { 0.15, 0.30, 0.50, 0.05 };
        em_match[3] = new double[] { 0.30, 0.15, 0.05, 0.50 };

        // Define emission for indels
        em_indel = 0.25;
    }

    public Pair<String, String> align(String query, String target) { return viterbi(query, target); }

    private int baseToIndex(char c) {
        switch (c) {
            case 'A': return 0;
            case 'T': return 1;
            case 'C': return 2;
            case 'G': return 3;
        }

        return 0;
    }

    private double emissionMatch(char a, char b) {
        int aIndex = baseToIndex(a);
        int bIndex = baseToIndex(b);

        return em_match[aIndex][bIndex];
    }

    private Pair<String, String> viterbi(String query, String target) {
        double[][] vm = new double[query.length() + 1][target.length() + 1];
        double[][] vi = new double[query.length() + 1][target.length() + 1];
        double[][] vd = new double[query.length() + 1][target.length() + 1];

        int[][] bm = new int[query.length() + 1][target.length() + 1];
        int[][] bi = new int[query.length() + 1][target.length() + 1];
        int[][] bd = new int[query.length() + 1][target.length() + 1];

        vm[0][0] = 1;
        for (int i = 1; i < query.length() + 1; i++) {
            for (int j = 1; j < target.length() + 1; j++) {
                if (i != 0 && j != 0) {
                    double s = emissionMatch(query.charAt(i - 1), target.charAt(j - 1));

                    double mm = tm[1][1] * vm[i - 1][j - 1];
                    double mi = tm[2][1] * vi[i - 1][j - 1];
                    double md = tm[3][1] * vd[i - 1][j - 1];

                    if (mm > mi && mm > md) {
                        bm[i][j] = 1;
                        vm[i][j] = s*mm;
                    } else if (mi > mm && mi > md) {
                        bm[i][j] = 2;
                        vm[i][j] = s*mi;
                    } else {
                        bm[i][j] = 3;
                        vm[i][j] = s*md;
                    }

                    double im = tm[1][2]*vm[i - 1][j];
                    double ii = tm[2][2]*vi[i - 1][j];

                    if (im > ii) {
                        bi[i][j] = 1;
                        vi[i][j] = em_indel * tm[1][2]*vm[i - 1][j];
                    } else {
                        bi[i][j] = 2;
                        vi[i][j] = em_indel * tm[2][2]*vm[i - 1][j];
                    }

                    double dm = tm[1][2]*vm[i][j - 1];
                    double dd = tm[3][3]*vd[i][j - 1];

                    if (dm > dd) {
                        bd[i][j] = 1;
                        vd[i][j] = em_indel * tm[1][3]*vm[i][j - 1];
                    } else {
                        bd[i][j] = 3;
                        vd[i][j] = em_indel * tm[3][3]*vm[i][j - 1];
                    }
                }
            }
        }

        StringBuilder qa = new StringBuilder();
        StringBuilder ta = new StringBuilder();

        int i = query.length();
        int j = target.length();
        int lastState = 0;

        switch (MoreMathUtils.whichMax(vm[i][j], vi[i][j], vd[i][j])) {
            case 0:
                qa.insert(0, query.charAt(i - 1));
                ta.insert(0, target.charAt(j - 1));
                lastState = bm[i][j];
                i -= 1;
                j -= 1;
                break;
            case 1:
                qa.insert(0, query.charAt(i - 1));
                ta.insert(0, "-");
                lastState = bi[i][j];
                i -= 1;
                break;
            case 2:
                qa.insert(0, "-");
                ta.insert(0, target.charAt(j - 1));
                lastState = bd[i][j];
                j -= 1;
                break;
        }

        while (i > 0 && j > 0) {
            switch (lastState) {
                case 1:
                    qa.insert(0, query.charAt(i - 1));
                    ta.insert(0, target.charAt(j - 1));
                    lastState = bm[i][j];
                    i -= 1;
                    j -= 1;
                    break;
                case 2:
                    qa.insert(0, query.charAt(i - 1));
                    ta.insert(0, "-");
                    lastState = bi[i][j];
                    i -= 1;
                    break;
                case 3:
                    qa.insert(0, "-");
                    ta.insert(0, target.charAt(j - 1));
                    lastState = bd[i][j];
                    j -= 1;
                    break;
            }
        }

        return new Pair<String, String>(qa.toString(), ta.toString());
    }

    private void printMatrix(String query, String target, double[][] m) {
        System.out.println("   " + Joiner.on(" ").join(Arrays.copyOfRange(target.split(""), 1, target.length() + 1)));

        for (int i = 0; i < m.length; i++) {
            List<Double> values = new ArrayList<Double>();
            for (int j = 0; j < m[0].length; j++) {
                values.add(m[i][j]);
            }

            if (i == 0) {
                System.out.println("^ " + Joiner.on(" ").join(values));
            } else {
                System.out.println(query.charAt(i - 1) + " " + Joiner.on(" ").join(values));
            }
        }
    }
}
