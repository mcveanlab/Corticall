package uk.ac.ox.well.cortexjdk.utils.alignment;

import org.testng.annotations.Test;
import uk.ac.ox.well.cortexjdk.utils.alignment.sw.SmithWaterman;

public class ConstrainedSmithWatermanTest {
    @Test(enabled = true)
    public void testAlignment() {
        String ref   = "TTTTCAGTACTAGTTCTAGTCATTATAAAACTGTCATG";
        String query = "GTCATTATAAAACTGTCATGTTTCATATATGTCATGGGGGCTAG";

        System.out.println(ref);
        System.out.println(query);

        /*
        String ref   = "ACACACTA";
        String query = "AGCACACA";

        ConstrainedSmithWaterman csw = new ConstrainedSmithWaterman(ref, query, 2, -1, -1);
        int[][] scoreTable = csw.getScoreTable();

        for (int[] aScoreTable : scoreTable) {
            for (int j = 0; j < aScoreTable.length; j++) {
                System.out.print(aScoreTable[j] + " ");
            }
            System.out.println();
        }

        String[] alignment = csw.getAlignment();
        System.out.println(alignment[0] + " " + alignment[1]);
        */

        SmithWaterman sw = new SmithWaterman(ref, query, 2, -1, -1);
        int[][] scoreTable2 = sw.getScoreTable();

        for (int[] aScoreTable2 : scoreTable2) {
            for (int j = 0; j < aScoreTable2.length; j++) {
                System.out.print(aScoreTable2[j] + " ");
            }
            System.out.println();
        }

        String[] alignment2 = sw.getAlignment();
        System.out.println(alignment2[0] + " " + alignment2[1]);
    }
}
