package uk.ac.ox.well.cortexjdk.utils.statistics;

import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.cortexjdk.utils.statistics.misc.StatisticsOnStream;

public class StatisticsOnStreamTest {
    StatisticsOnStream st = new StatisticsOnStream();

    public StatisticsOnStreamTest() {
        st.push(1);
        st.push(2);
        st.push(5);
        st.push(10);
        st.push(-5);
    }

    @Test
    public void getNumObservationsTest() {
        Assert.assertEquals(5, st.getNumObservations());
    }

    @Test
    public void getMeanTest() {
        Assert.assertEquals(2.6, st.getMean(), 0.01);
    }

    @Test
    public void getVarianceTest() {
        Assert.assertEquals(30.3, st.getVariance(), 0.01);
    }

    @Test
    public void getStandardDeviationTest() {
        Assert.assertEquals(5.504544, st.getStandardDeviation(), 0.01);
    }
}
