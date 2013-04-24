package uk.ac.ox.well.indiana.utils.alignment;

import org.testng.annotations.Test;

public class ConstrainedSmithWatermanTest {
    @Test
    public void testAlignment() {
        String ref = "ACACACTA";
        String query = "AGCACACA";

        System.out.println(ref);
        System.out.println(query);

        ConstrainedSmithWaterman csw = new ConstrainedSmithWaterman(ref, query, 2, -1, -1);
        int[][] scoreTable = csw.getScoreTable();

        for (int i = 0; i < scoreTable.length; i++) {
            for (int j = 0; j < scoreTable[i].length; j++) {
                System.out.print(scoreTable[i][j] + " ");
            }
            System.out.println();
        }

        String[] alignment = csw.getAlignment();
        System.out.println(alignment[0] + " " + alignment[1]);

        SmithWaterman sw = new SmithWaterman(ref, query, 2, -1, -1);
        int[][] scoreTable2 = sw.getScoreTable();

        for (int i = 0; i < scoreTable2.length; i++) {
            for (int j = 0; j < scoreTable2[i].length; j++) {
                System.out.print(scoreTable2[i][j] + " ");
            }
            System.out.println();
        }

        String[] alignment2 = sw.getAlignment();
        System.out.println(alignment2[0] + " " + alignment2[1]);
    }
}
