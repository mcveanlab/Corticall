package uk.ac.ox.well.indiana.utils.statistics.clustering;

import com.google.common.base.Joiner;
import picard.util.RExecutor;
import uk.ac.ox.well.indiana.utils.containers.DataFrame;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

public class HierarchicalClustering {
    private final String clusterGenesScript = "R/clustergenes.R";
    private DataFrame<String, String, Float> matrix;
    private boolean isSimilarityMatrix = true;
    private List<Set<String>> clusters = null;
    private int membershipLevel = -1;

    public void setMatrix(DataFrame<String, String, Float> m, boolean isSimilarityMatrix) {
        matrix = m;
        this.isSimilarityMatrix = isSimilarityMatrix;
    }

    public void setSimilarityMatrix(DataFrame<String, String, Float> m) { setMatrix(m, true); }
    public void setDissimilarityMatrix(DataFrame<String, String, Float> m) { setMatrix(m, false); }

    public void setMembershipLevel(int membershipLevel) { this.membershipLevel = membershipLevel; }
    public int getMembershipLevel() { return membershipLevel; }

    private void writeMatrix(File matrixTmpFile) {
        try {
            PrintStream out = new PrintStream(matrixTmpFile);

            out.println("\t" + Joiner.on("\t").join(matrix.getColNames()));

            for (String name1 : matrix.getRowNames()) {
                Collection<String> fields = new ArrayList<String>();
                fields.add(name1);

                for (String name2 : matrix.getColNames()) {
                    float value = matrix.get(name1, name2);
                    float distance = isSimilarityMatrix ? 1.0f - value : value;

                    fields.add(String.valueOf(distance));
                }

                out.println(Joiner.on("\t").join(fields));
            }

            out.close();
        } catch (FileNotFoundException e) {
            throw new IndianaException("Unable to write matrix file", e);
        }
    }

    private List<Set<String>> readClusters(File clustersFile) {
        List<Set<String>> clusters = new ArrayList<Set<String>>();

        TableReader tr = new TableReader(clustersFile);

        for (Map<String, String> te : tr) {
            String[] members = te.get("groupMembers").split(",");

            clusters.add(new HashSet<String>(Arrays.asList(members)));
        }

        return clusters;
    }

    public void cluster() {
        try {
            File matrixTmpFile = File.createTempFile("matrix", ".tmp");
            matrixTmpFile.deleteOnExit();
            writeMatrix(matrixTmpFile);

            //System.out.println(matrixTmpFile);

            File clustersFile = File.createTempFile("clusters", ".tmp");
            clustersFile.deleteOnExit();

            RExecutor.executeFromClasspath(clusterGenesScript, matrixTmpFile.getAbsolutePath(), clustersFile.getAbsolutePath(), String.valueOf(membershipLevel));

            clusters = readClusters(clustersFile);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public List<Set<String>> getClusters() { return clusters; }
}
