package uk.ac.ox.well.cortexjdk.utils.statistics.misc;

/**
 * Computes running statistics as each data point is added to the class.
 */
public class StatisticsOnStream {
    private int n = 0;

    double oldM, newM, oldS, newS;

    /**
     * Push a new data point onto the stream
     *
     * @param x  a data point to add
     */
    public void push(int x) {
        n++;

        if (n == 1) {
            oldM = newM = x;
            oldS = 0.0;
        } else {
            newM = oldM + (x - oldM)/n;
            newS = oldS + (x - oldM)*(x - newM);

            // set up for next iteration
            oldM = newM;
            oldS = newS;
        }
    }

    /**
     * Get the number of elements added to the stream in total.
     *
     * @return  number of observations
     */
    public int getNumObservations() { return n; }

    /**
     * Get the mean of the data added thus far.
     *
     * @return  mean of observations
     */
    public double getMean() { return (n > 0) ? newM : 0.0; }

    /**
     * Get the variance of the data added thus far.
     *
     * @return  variance of observations
     */
    public double getVariance() { return ( (n > 1) ? newS/(n - 1) : 0.0 ); }

    /**
     * Get the standard deviation of the data added thus far.
     *
     * @return  standard deviation of observations
     */
    public double getStandardDeviation() { return Math.sqrt( getVariance() ); }
}
