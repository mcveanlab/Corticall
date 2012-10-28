package uk.ac.ox.well.indiana.utils.statistics;

public class StatisticsOnStream {
    private int n = 0;

    double oldM, newM, oldS, newS;

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

    public int getNumObservations() { return n; }

    public double getMean() { return (n > 0) ? newM : 0.0; }

    public double getVariance() { return ( (n > 1) ? newS/(n - 1) : 0.0 ); }

    public double getStandardDeviation() { return Math.sqrt( getVariance() ); }
}
