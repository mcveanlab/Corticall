package uk.ac.ox.well.indiana.utils.math;

public class MoreMathUtils {
    private MoreMathUtils() {}

    public static double max(double... values) {
        double max = Double.NEGATIVE_INFINITY;

        for (double tmp : values) {
            max = max < tmp ? tmp : max;
        }

        return max;
    }

}
