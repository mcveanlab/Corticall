package uk.ac.ox.well.cortexjdk.utils.kmer;

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class CortexGraphTest {
    private class SimpleCortexRecord {
        public String kmer;
        public int[] coverage;
        public String[] edges;

        public SimpleCortexRecord(String kmer, int[] coverage, String[] edges) {
            this.kmer = kmer;
            this.coverage = coverage;
            this.edges = edges;
        }

        public boolean equals(CortexRecord cr) {
            if (kmer.equals(cr.getKmerAsString())) {
                if (coverage.length == cr.getCoverages().length) {
                    boolean coveragesAreEqual = true;
                    for (int i = 0; i < coverage.length; i++) {
                        coveragesAreEqual &= coverage[i] == cr.getCoverages()[i];
                    }

                    if (coveragesAreEqual) {
                        boolean edgesAreEqual = true;
                        for (int i = 0; i < edges.length; i++) {
                            edgesAreEqual &= edges[i].equals(cr.getEdgeAsStrings()[i]);
                        }

                        if (edgesAreEqual) {
                            return true;
                        }
                    }
                }
            }

            return false;
        }

        public String toString() {
            String value = kmer;

            for (int cov : coverage) {
                value += " " + cov;
            }

            for (String edge : edges) {
                value += " " + edge;
            }

            return value;
        }
    }

    @BeforeClass
    public void initialize() {
        recs = new ArrayList<>();

        recs.add(new SimpleCortexRecord("AAATAGGGCCACGATTTTTATTCAGAGCATA", new int[] { 1, 0 }, new String[] {"..g..C..", "........"}));
        recs.add(new SimpleCortexRecord("AACTACATGACCAGTACTCAGAGAGAAGCCC", new int[] { 0, 1 }, new String[] {"........", ".c..A..."}));
        recs.add(new SimpleCortexRecord("AATAGGGCCACGATTTTTATTCAGAGCATAC", new int[] { 1, 0 }, new String[] {"a.....G.", "........"}));
        recs.add(new SimpleCortexRecord("ACACTACGACTACAGCAACTACATGACCAGT", new int[] { 0, 1 }, new String[] {"........", "..g.A..."}));
        recs.add(new SimpleCortexRecord("ACAGCAACTACATGACCAGTACTCAGAGAGA", new int[] { 0, 1 }, new String[] {"........", "...tA..."}));
        recs.add(new SimpleCortexRecord("ACATGACCAGTACTCAGAGAGAAGCCCATAA", new int[] { 0, 1 }, new String[] {"........", "...t...T"}));
        recs.add(new SimpleCortexRecord("ACCAGTACTCAGAGAGAAGCCCATAATAGGC", new int[] { 0, 1 }, new String[] {"........", "..g...G."}));
        recs.add(new SimpleCortexRecord("ACGAAATAGGGCCACGATTTTTATTCAGAGC", new int[] { 1, 0 }, new String[] {"...tA...", "........"}));
        recs.add(new SimpleCortexRecord("ACGACACTACGACTACAGCAACTACATGACC", new int[] { 0, 1 }, new String[] {"........", ".c..A..."}));
        recs.add(new SimpleCortexRecord("ACGACTACAGCAACTACATGACCAGTACTCA", new int[] { 0, 1 }, new String[] {"........", "...t..G."}));
        recs.add(new SimpleCortexRecord("ACGATTTTTATTCAGAGCATACGATACAGAA", new int[] { 1, 0 }, new String[] {".c......", "........"}));
        recs.add(new SimpleCortexRecord("ACTACAGCAACTACATGACCAGTACTCAGAG", new int[] { 0, 1 }, new String[] {"........", "..g.A..."}));
        recs.add(new SimpleCortexRecord("ACTACATGACCAGTACTCAGAGAGAAGCCCA", new int[] { 0, 1 }, new String[] {"........", "a......T"}));
        recs.add(new SimpleCortexRecord("ACTACGACTACAGCAACTACATGACCAGTAC", new int[] { 0, 1 }, new String[] {"........", ".c.....T"}));
        recs.add(new SimpleCortexRecord("ACTATACGAAATAGGGCCACGATTTTTATTC", new int[] { 1, 0 }, new String[] {"....A...", "........"}));
        recs.add(new SimpleCortexRecord("ACTCAGAGAGAAGCCCATAATAGGCGCGGCC", new int[] { 0, 1 }, new String[] {"........", "...t.C.."}));
        recs.add(new SimpleCortexRecord("ACTGGGGGGCCACGACACTACGACTACAGCA", new int[] { 0, 1 }, new String[] {"........", "....A..."}));
        recs.add(new SimpleCortexRecord("AGCAACTACATGACCAGTACTCAGAGAGAAG", new int[] { 0, 1 }, new String[] {"........", ".c...C.."}));
        recs.add(new SimpleCortexRecord("AGGGCCACGATTTTTATTCAGAGCATACGAT", new int[] { 1, 0 }, new String[] {"...tA...", "........"}));
        recs.add(new SimpleCortexRecord("AGTACTCAGAGAGAAGCCCATAATAGGCGCG", new int[] { 0, 1 }, new String[] {"........", ".c....G."}));
        recs.add(new SimpleCortexRecord("AGTACTGGTCATGTAGTTGCTGTAGTCGTAG", new int[] { 0, 1 }, new String[] {"........", "..g....T"}));
        recs.add(new SimpleCortexRecord("AGTTGCTGTAGTCGTAGTGTCGTGGCCCCCC", new int[] { 0, 1 }, new String[] {"........", "...tA..."}));
        recs.add(new SimpleCortexRecord("ATACGAAATAGGGCCACGATTTTTATTCAGA", new int[] { 1, 0 }, new String[] {"...t..G.", "........"}));
        recs.add(new SimpleCortexRecord("ATAGGGCCACGATTTTTATTCAGAGCATACG", new int[] { 1, 0 }, new String[] {"a...A...", "........"}));
        recs.add(new SimpleCortexRecord("ATGACCAGTACTCAGAGAGAAGCCCATAATA", new int[] { 0, 1 }, new String[] {"........", ".c....G."}));
        recs.add(new SimpleCortexRecord("ATGCTCTGAATAAAAATCGTGGCCCTATTTC", new int[] { 1, 0 }, new String[] {"...t..G.", "........"}));
        recs.add(new SimpleCortexRecord("ATGGGCTTCTCTCTGAGTACTGGTCATGTAG", new int[] { 0, 1 }, new String[] {"........", "...t...T"}));
        recs.add(new SimpleCortexRecord("ATGTAGTTGCTGTAGTCGTAGTGTCGTGGCC", new int[] { 0, 1 }, new String[] {"........", ".c...C.."}));
        recs.add(new SimpleCortexRecord("ATTATGGGCTTCTCTCTGAGTACTGGTCATG", new int[] { 0, 1 }, new String[] {"........", "...t...T"}));
        recs.add(new SimpleCortexRecord("CAACTACATGACCAGTACTCAGAGAGAAGCC", new int[] { 0, 1 }, new String[] {"........", "..g..C.."}));
        recs.add(new SimpleCortexRecord("CACGACACTACGACTACAGCAACTACATGAC", new int[] { 0, 1 }, new String[] {"........", ".c...C.."}));
        recs.add(new SimpleCortexRecord("CACGATTTTTATTCAGAGCATACGATACAGA", new int[] { 1, 0 }, new String[] {".c..A...", "........"}));
        recs.add(new SimpleCortexRecord("CACTACGACTACAGCAACTACATGACCAGTA", new int[] { 0, 1 }, new String[] {"........", "a....C.."}));
        recs.add(new SimpleCortexRecord("CAGCAACTACATGACCAGTACTCAGAGAGAA", new int[] { 0, 1 }, new String[] {"........", "a.....G."}));
        recs.add(new SimpleCortexRecord("CAGTACTCAGAGAGAAGCCCATAATAGGCGC", new int[] { 0, 1 }, new String[] {"........", ".c....G."}));
        recs.add(new SimpleCortexRecord("CATGTAGTTGCTGTAGTCGTAGTGTCGTGGC", new int[] { 0, 1 }, new String[] {"........", "...t.C.."}));
        recs.add(new SimpleCortexRecord("CCACGACACTACGACTACAGCAACTACATGA", new int[] { 0, 1 }, new String[] {"........", "..g..C.."}));
        recs.add(new SimpleCortexRecord("CCACGATTTTTATTCAGAGCATACGATACAG", new int[] { 1, 0 }, new String[] {"..g.A...", "........"}));
        recs.add(new SimpleCortexRecord("CCAGTACTCAGAGAGAAGCCCATAATAGGCG", new int[] { 0, 1 }, new String[] {"........", "a....C.."}));
        recs.add(new SimpleCortexRecord("CCGCGCCTATTATGGGCTTCTCTCTGAGTAC", new int[] { 0, 1 }, new String[] {"........", "..g....T"}));
        recs.add(new SimpleCortexRecord("CCTATTATGGGCTTCTCTCTGAGTACTGGTC", new int[] { 0, 1 }, new String[] {"........", "..g.A..."}));
        recs.add(new SimpleCortexRecord("CGAAATAGGGCCACGATTTTTATTCAGAGCA", new int[] { 1, 0 }, new String[] {"a......T", "........"}));
        recs.add(new SimpleCortexRecord("CGACACTACGACTACAGCAACTACATGACCA", new int[] { 0, 1 }, new String[] {"........", "a.....G."}));
        recs.add(new SimpleCortexRecord("CGACTACAGCAACTACATGACCAGTACTCAG", new int[] { 0, 1 }, new String[] {"........", "a...A..."}));
        recs.add(new SimpleCortexRecord("CTACAGCAACTACATGACCAGTACTCAGAGA", new int[] { 0, 1 }, new String[] {"........", "a.....G."}));
        recs.add(new SimpleCortexRecord("CTATACGAAATAGGGCCACGATTTTTATTCA", new int[] { 1, 0 }, new String[] {"a.....G.", "........"}));
        recs.add(new SimpleCortexRecord("CTATTATGGGCTTCTCTCTGAGTACTGGTCA", new int[] { 0, 1 }, new String[] {"........", ".c.....T"}));
        recs.add(new SimpleCortexRecord("CTCAGAGAGAAGCCCATAATAGGCGCGGCCC", new int[] { 0, 1 }, new String[] {"........", "a......."}));
        recs.add(new SimpleCortexRecord("CTCTCTGAGTACTGGTCATGTAGTTGCTGTA", new int[] { 0, 1 }, new String[] {"........", "...t..G."}));
        recs.add(new SimpleCortexRecord("CTCTGAATAAAAATCGTGGCCCTATTTCGTA", new int[] { 1, 0 }, new String[] {"..g....T", "........"}));
        recs.add(new SimpleCortexRecord("CTGAATAAAAATCGTGGCCCTATTTCGTATA", new int[] { 1, 0 }, new String[] {"...t..G.", "........"}));
        recs.add(new SimpleCortexRecord("CTGGGGGGCCACGACACTACGACTACAGCAA", new int[] { 0, 1 }, new String[] {"........", "a....C.."}));
        recs.add(new SimpleCortexRecord("CTGGTCATGTAGTTGCTGTAGTCGTAGTGTC", new int[] { 0, 1 }, new String[] {"........", "a.....G."}));
        recs.add(new SimpleCortexRecord("GACTACAGCAACTACATGACCAGTACTCAGA", new int[] { 0, 1 }, new String[] {"........", ".c....G."}));
        recs.add(new SimpleCortexRecord("GAGTACTGGTCATGTAGTTGCTGTAGTCGTA", new int[] { 0, 1 }, new String[] {"........", "...t..G."}));
        recs.add(new SimpleCortexRecord("GCAACTACATGACCAGTACTCAGAGAGAAGC", new int[] { 0, 1 }, new String[] {"........", "a....C.."}));
        recs.add(new SimpleCortexRecord("GCCACGATTTTTATTCAGAGCATACGATACA", new int[] { 1, 0 }, new String[] {"..g...G.", "........"}));
        recs.add(new SimpleCortexRecord("GCCGCGCCTATTATGGGCTTCTCTCTGAGTA", new int[] { 0, 1 }, new String[] {"........", "..g..C.."}));
        recs.add(new SimpleCortexRecord("GGCCACGATTTTTATTCAGAGCATACGATAC", new int[] { 1, 0 }, new String[] {"..g.A...", "........"}));
        recs.add(new SimpleCortexRecord("GGGCCACGACACTACGACTACAGCAACTACA", new int[] { 0, 1 }, new String[] {"........", "..g....T"}));
        recs.add(new SimpleCortexRecord("GGGCCACGATTTTTATTCAGAGCATACGATA", new int[] { 1, 0 }, new String[] {"a....C..", "........"}));
        recs.add(new SimpleCortexRecord("GGGGCCACGACACTACGACTACAGCAACTAC", new int[] { 0, 1 }, new String[] {"........", "..g.A..."}));
        recs.add(new SimpleCortexRecord("GGGGGCCACGACACTACGACTACAGCAACTA", new int[] { 0, 1 }, new String[] {"........", "..g..C.."}));
        recs.add(new SimpleCortexRecord("GTTGCTGTAGTCGTAGTGTCGTGGCCCCCCA", new int[] { 0, 1 }, new String[] {"........", "a.....G."}));
        recs.add(new SimpleCortexRecord("TACATGACCAGTACTCAGAGAGAAGCCCATA", new int[] { 0, 1 }, new String[] {"........", ".c..A..."}));
        recs.add(new SimpleCortexRecord("TAGGGCCACGATTTTTATTCAGAGCATACGA", new int[] { 1, 0 }, new String[] {"a......T", "........"}));
    }

    @Test
    public void getShortSampleNames() {
        CortexGraph cg = new CortexGraph("testdata/two_short_contigs.ctx");

        Assert.assertEquals("one", cg.getColor(0).getSampleName());
        Assert.assertEquals("two", cg.getColor(1).getSampleName());
    }

    @Test
    public void numRecordsTest() {
        CortexGraph cg = new CortexGraph("testdata/two_short_contigs.ctx");

        Assert.assertEquals(66, cg.getNumRecords());
    }

    @Test
    public void getLeftAndRightEdges() {
        CortexGraph cg = new CortexGraph("testdata/two_short_contigs.ctx");

        CortexRecord cr = cg.next();

        for (int color = 0; color < cg.getNumColors(); color++) {
            byte[] edges = new String(cr.getEdgesAsBytes(color)).toUpperCase().getBytes();

            Set<Byte> leftEdges = new HashSet<>();
            Set<Byte> rightEdges = new HashSet<>();

            for (int i = 0; i < 4; i++) {
                if (edges[i] != '.') {
                    leftEdges.add(edges[i]);
                }

                if (edges[i+4] != '.') {
                    rightEdges.add(edges[i+4]);
                }
            }

            Set<Byte> betterLeftEdges = new HashSet<>(cr.getInEdgesAsBytes(color, false));
            Set<Byte> betterRightEdges = new HashSet<>(cr.getOutEdgesAsBytes(color, false));

            Assert.assertEquals(leftEdges, betterLeftEdges);
            Assert.assertEquals(rightEdges, betterRightEdges);
        }
    }

    private List<SimpleCortexRecord> recs;

    @Test
    public void recordsAreCorrect() {
        CortexGraph cg = new CortexGraph("testdata/two_short_contigs.ctx");

        int index = 0;
        for (CortexRecord cr : cg) {
            SimpleCortexRecord scr = recs.get(index);

            Assert.assertEquals(scr.equals(cr), true, "Cortex record says '" + cr + "' but test record says '" + scr + "'");

            index++;
        }
    }

    @Test
    public void constructRecordsFromString() {
        CortexGraph cg = new CortexGraph("testdata/two_short_contigs.ctx");

        for (CortexRecord cr : cg) {
            CortexRecord nr = new CortexRecord(cr.toString());

            Assert.assertEquals(nr, cr);
        }
    }

    @Test
    public void constructRecords() {
        CortexGraph cg = new CortexGraph("testdata/two_short_contigs.ctx");

        for (CortexRecord cr : cg) {
            String sk = cr.getKmerAsString();
            List<Integer> coverageList = new ArrayList<>();
            List<Set<String>> inEdgesList = new ArrayList<>();
            List<Set<String>> outEdgesList = new ArrayList<>();

            for (int c = 0; c < cr.getNumColors(); c++) {
                coverageList.add(cr.getCoverage(c));
                inEdgesList.add(new HashSet<>(cr.getInEdgesAsStrings(c, false)));
                outEdgesList.add(new HashSet<>(cr.getOutEdgesAsStrings(c, false)));
            }

            CortexRecord nr = new CortexRecord(sk, coverageList, inEdgesList, outEdgesList);

            Assert.assertEquals(nr, cr);
        }
    }

    @Test
    public void constructRcRecords() {
        CortexGraph cg = new CortexGraph("testdata/two_short_contigs.ctx");

        for (CortexRecord cr : cg) {
            String sk = SequenceUtils.reverseComplement(cr.getKmerAsString());
            List<Integer> coverageList = new ArrayList<>();
            List<Set<String>> inEdgesList = new ArrayList<>();
            List<Set<String>> outEdgesList = new ArrayList<>();

            for (int c = 0; c < cr.getNumColors(); c++) {
                coverageList.add(cr.getCoverage(c));
                inEdgesList.add(new HashSet<>(cr.getOutEdgesAsStrings(c, true)));
                outEdgesList.add(new HashSet<>(cr.getInEdgesAsStrings(c, true)));
            }

            CortexRecord nr = new CortexRecord(sk, coverageList, inEdgesList, outEdgesList);

            Assert.assertEquals(nr, cr);
        }
    }

    @Test
    public void testGetRecord() {
        CortexGraph cg = new CortexGraph("testdata/two_short_contigs.ctx");

        for (int i = 10; i >= 0; i--) {
            CortexRecord cr = cg.getRecord(i);
            SimpleCortexRecord scr = recs.get(i);

            Assert.assertEquals(scr.equals(cr), true, "Cortex record says '" + cr + "' but test record says '" + scr + "'");
        }
    }

    @Test
    public void testEncodeBinaryKmer() {
        CortexGraph cg = new CortexGraph("testdata/two_short_contigs.ctx");

        for (int i = 10; i >= 0; i--) {
            CortexRecord cr = cg.getRecord(i);

            long[] originalBinaryKmer = cr.getBinaryKmer();
            long[] encodedBinaryKmer = CortexRecord.encodeBinaryKmer(cr.getKmerAsBytes());

            Assert.assertEquals(encodedBinaryKmer, originalBinaryKmer);
            Assert.assertEquals(CortexRecord.decodeBinaryKmer(encodedBinaryKmer, cr.getKmerSize(), cr.getKmerBits()), cr.getKmerAsBytes());
        }
    }

    /*
    @Test
    public void getBigKmerTest() {
        CortexGraph cg = new CortexGraph("testdata/Pf3D7_01_v3.k95.ctx");

        for (int i = 0; i < 100; i++) {
            CortexRecord cr = cg.getRecord(i);

            byte[] kmer = cr.getKmerAsBytes();
            long[] oldbk = cr.getBinaryKmer();
            long[] newbk = CortexRecord.encodeBinaryKmer(kmer);

            Assert.assertEquals(oldbk, newbk);
            Assert.assertEquals(kmer, CortexRecord.decodeBinaryKmer(newbk, cr.getKmerSize(), cr.getKmerBits()));
        }
    }
    */

    /*
    @Test(expectedExceptions = CortexJDKException.class)
    public void testUnsortedFindRecordThrowsException() {
        CortexGraph cg = new CortexGraph("testdata/two_short_contigs.ctx");
        String targetKmer = "ACCAAATCATCTGTTGGAAATTTATTACAAA";

        cg.findRecord(targetKmer);
    }
    */

    @Test
    public void testSortedFindRecord() {
        CortexGraph cg = new CortexGraph("testdata/two_short_contigs.ctx");

        for (SimpleCortexRecord scr : recs) {
            CortexRecord cr = cg.findRecord(scr.kmer);

            Assert.assertNotNull(cr);
            Assert.assertEquals(scr.equals(cr), true);
        }
    }

    @Test
    public void testFindNonExistentRecord() {
        CortexGraph cg = new CortexGraph("testdata/two_short_contigs.ctx");

        String nonExistentRecord = "NTTTTGGGGTATTTGCAGTATTTGGAATAAA";

        CortexRecord cr = cg.findRecord(nonExistentRecord);

        Assert.assertNull(cr);
    }
}
