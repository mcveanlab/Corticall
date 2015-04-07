package uk.ac.ox.well.indiana.utils.alignment.pairwise;

import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.indiana.utils.math.MoreMathUtils;

public class LocalAligner {
    //private String[] states = new String[] { "M", "I", "D" };
    //private String[] emissions = new String[] { "+1,+1", "+1,+0", "+0,+1" };

    private double delta;
    private double epsilon;

    private double[][] tm;
    private double[][] em;

    public LocalAligner() {
        initialize(0.2, 0.2);
    }

    public LocalAligner(double delta, double epsilon) {
        initialize(delta, epsilon);
    }

    private void initialize(double delta, double epsilon) {
        this.delta = delta;
        this.epsilon = epsilon;

        // Define transition matrix
        tm = new double[3][3];
        tm[0] = new double[] { 1 - 2*delta, delta,   delta   };
        tm[1] = new double[] { 1 - epsilon, epsilon, 0       };
        tm[2] = new double[] { 1 - epsilon, 0,       epsilon };

        // Define emission matrix
        em = new double[3][3];
        em[0] = new double[] { 1, 0, 0 };
        em[1] = new double[] { 0, 1, 0 };
        em[2] = new double[] { 0, 0, 1 };
    }

    public Pair<String, String> align(String query, String target) {
        return viterbi(query, target, tm, em);
    }

    private Pair<String, String> viterbi(String query, String target, double[][] tm, double[][] em) {
        double[][] vm = new double[query.length()][target.length()];
        double[][] vi = new double[query.length()][target.length()];
        double[][] vd = new double[query.length()][target.length()];

        vm[0][0] = 1;
        for (int i = 0; i < query.length(); i++) {
            for (int j = 0; j < target.length(); j++) {
                if (i != 0 && j != 0) {
                    vm[i][j] = em[0][0] * MoreMathUtils.max(
                            (1 - 2 * delta) * vm[i - 1][j - 1],
                            (1 - epsilon) * vi[i - 1][j - 1],
                            (1 - epsilon) * vd[i - 1][j - 1]
                    );

                    vi[i][j] = em[1][1] * MoreMathUtils.max(
                            (delta)*vm[i - 1][j],
                            (epsilon)*vi[i - 1][j]
                    );

                    vd[i][j] = em[2][2] * MoreMathUtils.max(
                            (delta)*vm[i][j - 1],
                            (epsilon)*vd[i][j - 1]
                    );
                }
            }
        }

        StringBuilder qa = new StringBuilder();
        StringBuilder ta = new StringBuilder();

        int i = query.length() - 1;
        int j = target.length() - 1;

        while (i >= 0 && j >= 0) {
            switch (MoreMathUtils.whichMax(vm[i][j], vi[i][j], vd[i][j])) {
                case 0:
                    qa.insert(0, query.charAt(i));
                    ta.insert(0, target.charAt(j));
                    i -= 1;
                    j -= 1;
                    break;
                case 1:
                    qa.insert(0, query.charAt(i));
                    ta.insert(0, "-");
                    i -= 1;
                    break;
                case 2:
                    qa.insert(0, "-");
                    ta.insert(0, target.charAt(j));
                    j -= 1;
                    break;
            }
        }

        return new Pair<String, String>(qa.toString(), ta.toString());
    }
}
