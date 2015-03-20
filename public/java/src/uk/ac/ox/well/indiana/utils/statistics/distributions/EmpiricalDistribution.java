package uk.ac.ox.well.indiana.utils.statistics.distributions;

import uk.ac.ox.well.indiana.utils.math.MoreMathUtils;

import java.util.Random;

public class EmpiricalDistribution {
    private double[] rates;
    private double[] cdf;
    private final double tol = 1e-8;

    private Random rng;

    public EmpiricalDistribution(int[] dist) {
        initialize(dist, new Random(System.currentTimeMillis()));
    }

    public EmpiricalDistribution(int[] dist, long seed) {
        initialize(dist, new Random(seed));
    }

    public EmpiricalDistribution(int[] dist, Random rng) {
        initialize(dist, rng);
    }

    private void initialize(int[] dist, Random rng) {
        this.rng = rng;

        this.rates = new double[dist.length];

        long sum = 0;
        for (int val : dist) {
            sum += val;
        }

        for (int i = 0; i < dist.length; i++) {
            rates[i] = (double) dist[i] / (double) sum;
        }

        double cdfSum = 0.0;
        cdf = new double[rates.length];
        for (int i = 0; i < rates.length; i++) {
            cdfSum += rates[i];
            cdf[i] = cdfSum;
        }
    }

    public double[] getRates() {
        return rates;
    }

    private int findBin(double cdfValue) {
        int lowerBin = 0;
        int higherBin = cdf.length - 1;
        int middleBin = lowerBin + Math.round(((float) (higherBin - lowerBin)) / 2.0f);

        while (lowerBin != middleBin && middleBin != higherBin) {
            double lowerCdfValue = cdf[lowerBin];
            double middleCdfValue = cdf[middleBin];
            double higherCdfValue = cdf[higherBin];

            if (MoreMathUtils.equals(cdfValue, lowerCdfValue, tol)) {
                return lowerBin;
            } else if (MoreMathUtils.equals(cdfValue, middleCdfValue, tol)) {
                return middleBin;
            } else if (MoreMathUtils.equals(cdfValue, higherCdfValue, tol)) {
                return higherBin;
            } else if (lowerCdfValue < cdfValue && cdfValue < middleCdfValue) {
                higherBin = middleBin;
                middleBin = lowerBin + Math.round(((float) (higherBin - lowerBin)) / 2.0f);
            } else if (middleCdfValue < cdfValue && cdfValue < higherCdfValue) {
                lowerBin = middleBin;
                middleBin = lowerBin + Math.round(((float) (higherBin - lowerBin)) / 2.0f);
            } else if (cdfValue < lowerCdfValue) {
                return lowerBin;
            } else if (cdfValue > higherCdfValue) {
                return higherBin;
            } else {
                return middleBin;
            }
        }

        return middleBin;
    }

    public int draw() {
        return findBin(rng.nextDouble());
    }
}
