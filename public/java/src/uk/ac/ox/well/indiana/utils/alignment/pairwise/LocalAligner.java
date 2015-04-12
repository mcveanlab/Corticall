package uk.ac.ox.well.indiana.utils.alignment.pairwise;

import com.google.common.base.Joiner;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.math.MoreMathUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class LocalAligner {
    private double[][] tm;       // transition matrix
    private double[][] em_match; // emission matrix (match/mismatches)
    private double em_indel;     // emission matrix (indels)

    private double[][][] v;
    private int[][][] b;

    public LocalAligner() { initialize(0.2, 0.1, 0.1, 1e-4); }

    public LocalAligner(double delta, double epsilon, double tau, double eta) { initialize(delta, epsilon, tau, eta); }

    private void initialize(double delta, double epsilon, double tau, double eta) {
        // Define transition matrix
        tm = new double[13][13];

        //                      b  rx1,     fl,  ry1,     ab,  M,                 X,       Y,       ae,      rx2,     fr,  ry2,     e
        tm[0]  = new double[] { 0, 1 - eta, eta, 0,       0,   0,                 0,       0,       0,       0,       0,   0,       0   }; // b
        tm[1]  = new double[] { 0, 1 - eta, eta, 0,       0,   0,                 0,       0,       0,       0,       0,   0,       0   }; // rx1
        tm[2]  = new double[] { 0, 0,       0,   1 - eta, eta, 0,                 0,       0,       0,       0,       0,   0,       0   }; // fl
        tm[3]  = new double[] { 0, 0,       0,   1 - eta, eta, 0,                 0,       0,       0,       0,       0,   0,       0   }; // ry1

        tm[4]  = new double[] { 0, 0,       0,   0,       0,   1 - 2*delta - tau, delta,   delta,   tau,     0,       0,   0,       0   }; // ab
        tm[5]  = new double[] { 0, 0,       0,   0,       0,   1 - 2*delta - tau, delta,   delta,   tau,     0,       0,   0,       0   }; // M
        tm[6]  = new double[] { 0, 0,       0,   0,       0,   1 - epsilon - tau, epsilon, 0,       tau,     0,       0,   0,       0   }; // X
        tm[7]  = new double[] { 0, 0,       0,   0,       0,   1 - epsilon - tau, 0,       epsilon, tau,     0,       0,   0,       0   }; // Y
        tm[8]  = new double[] { 0, 0,       0,   0,       0,   0,                 0,       0,       0,       1 - eta, eta, 0,       0   }; // ae

        tm[9]  = new double[] { 0, 0,       0,   0,       0,   0,                 0,       0,       0,       1 - eta, eta, 0,       0   }; // rx2
        tm[10] = new double[] { 0, 0,       0,   0,       0,   0,                 0,       0,       0,       0,       0,   1 - eta, eta }; // fr
        tm[11] = new double[] { 0, 0,       0,   0,       0,   0,                 0,       0,       0,       0,       0,   1 - eta, eta }; // ry2
        tm[12] = new double[] { 0, 0,       0,   0,       0,   0,                 0,       0,       0,       0,       0,   0,       0   }; // e

        // Define emission for matches
        em_match = new double[4][4];
        em_match[0] = new double[] { 0.85, 0.05, 0.05, 0.05 };
        em_match[1] = new double[] { 0.05, 0.85, 0.05, 0.05 };
        em_match[2] = new double[] { 0.05, 0.05, 0.85, 0.05 };
        em_match[3] = new double[] { 0.05, 0.05, 0.05, 0.85 };

        // Define emission for indels
        em_indel = 0.25;
    }

    public Pair<String, String> align(String query, String target) {
        viterbi(query, target);
        return traceback(query, target);
    }

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

    private void viterbi(String query, String target) {
        v = new double[tm.length][query.length() + 1][target.length() + 1];
        b = new int[tm.length][query.length() + 1][target.length() + 1];

        for (int s = 0; s < tm.length; s++) {
            for (int i = 0; i < query.length() + 1; i++) {
                for (int j = 0; j < target.length() + 1; j++) {
                    v[s][i][j] = 0.0;
                }
            }
        }

        v[0][0][0] = 1;

        for (int i = 1; i < query.length() + 1; i++) {
            for (int j = 1; j < target.length() + 1; j++) {
                if (i != 0 && j != 0) {
                    double pab  = emissionMatch(query.charAt(i - 1), target.charAt(j - 1));

                    for (int sa = 0; sa < tm.length; sa++) {
                        double[] sb = new double[tm.length];

                        //sb[sa] =
                    }

                    /*

                    double s = Math.log10(pab);

                    double mm = s + Math.log10(tm[1][1]) + vm[i - 1][j - 1];
                    double mi = s + Math.log10(tm[2][1]) + vi[i - 1][j - 1];
                    double md = s + Math.log10(tm[3][1]) + vd[i - 1][j - 1];

                    if (mm > mi && mm > md) {
                        bm[i][j] = 1;
                        vm[i][j] = mm;
                    } else if (mi > mm && mi > md) {
                        bm[i][j] = 2;
                        vm[i][j] = mi;
                    } else {
                        bm[i][j] = 3;
                        vm[i][j] = md;
                    }

                    double im = Math.log10(em_indel) + Math.log10(tm[1][1]) + vm[i - 1][j];
                    double ii = Math.log10(em_indel) + Math.log10(tm[1][2]) + vi[i - 1][j];

                    if (im > ii) {
                        bi[i][j] = 1;
                        vi[i][j] = im;
                    } else {
                        bi[i][j] = 2;
                        vi[i][j] = ii;
                    }

                    double dm = Math.log10(em_indel) + Math.log10(tm[1][1]) + vm[i][j - 1];
                    double dd = Math.log10(em_indel) + Math.log10(tm[1][3]) + vd[i][j - 1];

                    if (dm > dd) {
                        bd[i][j] = 1;
                        vd[i][j] = dm;
                    } else {
                        bd[i][j] = 3;
                        vd[i][j] = dd;
                    }
                    */
                }
            }
        }
    }

    private Pair<String, String> traceback(String query, String target) {
        StringBuilder qa = new StringBuilder();
        StringBuilder ta = new StringBuilder();

        int i = query.length();
        int j = target.length();
        int lastState = 0;

        /*
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
            default:
                throw new IndianaException("Saw weird max");
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
                default:
                    throw new IndianaException("Saw weird state: " + lastState);
            }
        }
        */

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
