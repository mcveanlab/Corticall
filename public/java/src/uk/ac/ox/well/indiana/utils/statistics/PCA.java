package uk.ac.ox.well.indiana.utils.statistics;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import uk.ac.ox.well.indiana.utils.containers.DataFrame;

import java.util.HashMap;

public class PCA<R, C> {
    private SingularValueDecomposition svd;

    private HashMap<R, Integer> rowMapping = new HashMap<R, Integer>();
    private HashMap<C, Integer> colMapping = new HashMap<C, Integer>();

    public PCA(DataFrame<R, C, Float> d) {
        DoubleMatrix2D m = DoubleFactory2D.dense.make(d.getNumRows(), d.getNumCols());

        int r = 0;
        for (R rowName : d.getRowNames()) {
            int c = 0;

            for (C colName : d.getColNames()) {
                float datum = d.get(rowName, colName);

                m.set(r, c, datum);

                colMapping.put(colName, c);

                c++;
            }

            rowMapping.put(rowName, r);

            r++;
        }

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

        svd = new SingularValueDecomposition(m);
    }

    public SingularValueDecomposition getSVD() {
        return svd;
    }

    public DoubleMatrix1D getEigenvector(C colName) {
        return svd.getV().viewRow(colMapping.get(colName));
    }
}
