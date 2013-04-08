package uk.ac.ox.well.indiana.utils.statistics;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import uk.ac.ox.well.indiana.utils.containers.DataFrame;

import java.util.HashMap;

/**
 * Perform a principal components analysis on a dataset.
 *
 * @param <R>  The data type for the rows
 * @param <C>  The data type for the columns
 */
public class PCA<R, C> {
    private SingularValueDecomposition svd;

    private HashMap<R, Integer> rowMapping = new HashMap<R, Integer>();
    private HashMap<C, Integer> colMapping = new HashMap<C, Integer>();

    /**
     * Compute a PCA on the provided data frame.
     *
     * @param d  the data frame upon which the PCA should be performed
     */
    public PCA(DataFrame<R, C, Float> d) {
        computePCA(d, true, true);
    }

    /**
     * Compute a PCA on the provided data frame, optionally centering and scaling the data.
     *
     * @param d  the data frame upon which the PCA should be performed
     * @param center  true if the data should be centered, false otherwise
     * @param scale  true if the data should be scaled, false otherwise
     */
    public PCA(DataFrame<R, C, Float> d, boolean center, boolean scale) {
        computePCA(d, center, scale);
    }

    /**
     * Compute the actual PCA.
     *
     * @param d  the data frame upon which the PCA should be performed
     * @param center  true if the data should be centered, false otherwise
     * @param scale  true if the data should be scaled, false otherwise
     */
    private void computePCA(DataFrame<R, C, Float> d, boolean center, boolean scale) {
        DoubleMatrix2D m = DoubleFactory2D.sparse.make(d.getNumRows(), d.getNumCols());

        int r = 0;
        for (R rowName : d.getRowNames()) {
            int c = 0;

            for (C colName : d.getColNames()) {
                if (d.hasValue(rowName, colName)) {
                    float datum = d.get(rowName, colName);

                    m.set(r, c, datum);
                }

                colMapping.put(colName, c);

                c++;
            }

            rowMapping.put(rowName, r);

            r++;
        }

        if (center) {
            for (int c1 = 0; c1 < m.columns(); c1++) {
                double sum = 0.0;
                for (int r1 = 0; r1 < m.rows(); r1++) {
                    sum += m.get(r1, c1);
                }

                double mean = sum / ((double) m.rows());

                for (int r1 = 0; r1 < m.rows(); r1++) {
                    m.set(r1, c1, m.get(r1, c1) - mean);
                }
            }
        }

        if (scale) {
            for (int c1 = 0; c1 < m.columns(); c1++) {
                double mean = 0.0;

                if (center) {
                    double sum = 0.0;
                    for (int r1 = 0; r1 < m.rows(); r1++) {
                        sum += m.get(r1, c1);
                    }

                    mean = sum / ((double) m.rows());
                }

                double summedSquares = 0.0;

                for (int r1 = 0; r1 < m.rows(); r1++) {
                    summedSquares += Math.pow(m.get(r1, c1) - mean, 2.0);
                }

                double sd = Math.sqrt( summedSquares / ((double) m.rows()) );

                for (int r1 = 0; r1 < m.rows(); r1++) {
                    m.set(r1, c1, m.get(r1, c1)/sd);
                }
            }
        }

        svd = new SingularValueDecomposition(m);
    }

    /**
     * Return an eigenvector associated with the specified column.
     *
     * @param colName  the column name for which the eigenvector should be returned
     * @return  the eigenvector
     */
    public DoubleMatrix1D getEigenvector(C colName) {
        return svd.getV().viewRow(colMapping.get(colName));
    }

    /**
     * Return the rotation matrix of the PCA.
     *
     * @return  the rotation matrix
     */
    public DataFrame<C, String, Float> getRotationMatrix() {
        DataFrame<C, String, Float> d = new DataFrame<C, String, Float>(0.0f);

        for (C colName : colMapping.keySet()) {
            DoubleMatrix1D dm = svd.getV().viewRow(colMapping.get(colName));

            for (int pc = 0; pc < dm.size(); pc++) {
                float value = (float) dm.get(pc);

                String pcName = "PC" + (pc+1);
                d.set(colName, pcName, value);
            }
        }

        return d;
    }

    /**
     * Return the standard deviations of the PCA.
     *
     * @return  the standard deviations
     */
    public double[] getStandardDeviations() {
        double[] sds = svd.getSingularValues();

        for (int i = 0; i < sds.length; i++) {
            sds[i] /= Math.sqrt(Math.max(1.0, rowMapping.size() - 1));
        }

        return sds;
    }
}
