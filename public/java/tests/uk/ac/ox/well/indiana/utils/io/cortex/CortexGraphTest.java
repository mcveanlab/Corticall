package uk.ac.ox.well.indiana.utils.io.cortex;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;

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
            if (kmer.equals(cr.getKmerString())) {
                if (coverage.length == cr.getCoverages().length) {
                    boolean coveragesAreEqual = true;
                    for (int i = 0; i < coverage.length; i++) {
                        coveragesAreEqual &= coverage[i] == cr.getCoverages()[i];
                    }

                    if (coveragesAreEqual) {
                        boolean edgesAreEqual = true;
                        for (int i = 0; i < edges.length; i++) {
                            edgesAreEqual &= edges[i].equals(cr.getEdges()[i]);
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

    @Test
    public void getLongSampleNames() {
        CortexGraph cg = new CortexGraph("testdata/test_gene_for_sn_reconstruction.ctx");

        Assert.assertEquals("test_sample_with_long_name", cg.getColor(0).getSampleName());
    }

    @Test
    public void getShortSampleNames() {
        CortexGraph cg = new CortexGraph("testdata/smallgraph.ctx");

        Assert.assertEquals("ref", cg.getColor(0).getSampleName());
        Assert.assertEquals("ref_ss", cg.getColor(1).getSampleName());
    }

    @Test
    public void numRecordsTest() {
        CortexGraph cg = new CortexGraph("testdata/smallgraph.ctx");

        Assert.assertEquals(85, cg.getNumRecords());
    }

    @Test
    public void recordsAreCorrect() {
        ArrayList<SimpleCortexRecord> recs = new ArrayList<SimpleCortexRecord>();

        recs.add(new SimpleCortexRecord("GTATTTGCAGTATTTGGAATAAATTTCCAAC", new int[] {1,0}, new String[] {"..g.A...", "........"}));
        recs.add(new SimpleCortexRecord("GGAAATTTATTCCAAATACTGCAAATACCCC", new int[] {1,0}, new String[] {"...tA...", "........"}));
        recs.add(new SimpleCortexRecord("GTATTTGCAGTATTTGTAATAAATTTCCAAC", new int[] {0,1}, new String[] {"........", "..g.A..."}));
        recs.add(new SimpleCortexRecord("AAATCATCTGTTGGAAATTTATTACAAATAC", new int[] {0,1}, new String[] {"........", ".c.....T"}));
        recs.add(new SimpleCortexRecord("GGTATTTGCAGTATTTGTAATAAATTTCCAA", new int[] {0,1}, new String[] {"........", "..g..C.."}));
        recs.add(new SimpleCortexRecord("AAAAAACCAAATCATCTGTTGGAAATTTATT", new int[] {1,1}, new String[] {"a....C..", "a...A..."}));
        recs.add(new SimpleCortexRecord("AAAATCGTTTTGGGGTATTTGCAGTATTTGT", new int[] {0,1}, new String[] {"........", "..g.A..."}));
        recs.add(new SimpleCortexRecord("GCAAATACCCCAAAACGATTTTCTATAGCTA", new int[] {1,1}, new String[] {"...t...T", "...t...T"}));
        recs.add(new SimpleCortexRecord("ATCATCTGTTGGAAATTTATTACAAATACTG", new int[] {0,1}, new String[] {"........", "a....C.."}));
        recs.add(new SimpleCortexRecord("AAACCAAATCATCTGTTGGAAATTTATTACA", new int[] {0,1}, new String[] {"........", "a...A..."}));
        recs.add(new SimpleCortexRecord("AACCAAATCATCTGTTGGAAATTTATTACAA", new int[] {0,1}, new String[] {"........", "a...A..."}));
        recs.add(new SimpleCortexRecord("CGTTTTGGGGTATTTGCAGTATTTGGAATAA", new int[] {1,0}, new String[] {"...tA...", "........"}));
        recs.add(new SimpleCortexRecord("CTGCAAATACCCCAAAACGATTTTCTATAGC", new int[] {1,1}, new String[] {"a......T", "a......T"}));
        recs.add(new SimpleCortexRecord("GAAATTTATTCCAAATACTGCAAATACCCCA", new int[] {1,0}, new String[] {"..g.A...", "........"}));
        recs.add(new SimpleCortexRecord("AAATACCCCAAAACGATTTTCTATAGCTATG", new int[] {1,1}, new String[] {".c.....T", ".c.....T"}));
        recs.add(new SimpleCortexRecord("CAAATCATCTGTTGGAAATTTATTCCAAATA", new int[] {1,0}, new String[] {".c...C..", "........"}));
        recs.add(new SimpleCortexRecord("ATTTATTACAAATACTGCAAATACCCCAAAA", new int[] {0,1}, new String[] {"........", "a....C.."}));
        recs.add(new SimpleCortexRecord("AAAAAAACCAAATCATCTGTTGGAAATTTAT", new int[] {1,1}, new String[] {".......T", ".......T"}));
        recs.add(new SimpleCortexRecord("AAACCAAATCATCTGTTGGAAATTTATTCCA", new int[] {1,0}, new String[] {"a...A...", "........"}));
        recs.add(new SimpleCortexRecord("AATACCCCAAAACGATTTTCTATAGCTATGT", new int[] {1,1}, new String[] {"a...A...", "a...A..."}));
        recs.add(new SimpleCortexRecord("ATTTGCAGTATTTGGAATAAATTTCCAACAG", new int[] {1,0}, new String[] {"...tA...", "........"}));
        recs.add(new SimpleCortexRecord("ATCGTTTTGGGGTATTTGCAGTATTTGGAAT", new int[] {1,0}, new String[] {"a...A...", "........"}));
        recs.add(new SimpleCortexRecord("AAAATCGTTTTGGGGTATTTGCAGTATTTGG", new int[] {1,0}, new String[] {"..g.A...", "........"}));
        recs.add(new SimpleCortexRecord("GCAGTATTTGGAATAAATTTCCAACAGATGA", new int[] {1,0}, new String[] {"...t...T", "........"}));
        recs.add(new SimpleCortexRecord("ATACCCCAAAACGATTTTCTATAGCTATGTA", new int[] {1,1}, new String[] {"a.....G.", "a.....G."}));
        recs.add(new SimpleCortexRecord("GTTTTGGGGTATTTGCAGTATTTGGAATAAA", new int[] {1,0}, new String[] {".c.....T", "........"}));
        recs.add(new SimpleCortexRecord("AAATCGTTTTGGGGTATTTGCAGTATTTGTA", new int[] {0,1}, new String[] {"........", "a...A..."}));
        recs.add(new SimpleCortexRecord("ATTTGTAATAAATTTCCAACAGATGATTTGG", new int[] {0,1}, new String[] {"........", "...t...T"}));
        recs.add(new SimpleCortexRecord("ATTTATTCCAAATACTGCAAATACCCCAAAA", new int[] {1,0}, new String[] {"a....C..", "........"}));
        recs.add(new SimpleCortexRecord("TCTGTTGGAAATTTATTACAAATACTGCAAA", new int[] {0,1}, new String[] {"........", "a......T"}));
        recs.add(new SimpleCortexRecord("AAATCGTTTTGGGGTATTTGCAGTATTTGGA", new int[] {1,0}, new String[] {"a...A...", "........"}));
        recs.add(new SimpleCortexRecord("CAAATACTGCAAATACCCCAAAACGATTTTC", new int[] {1,1}, new String[] {".c.....T", "a......T"}));
        recs.add(new SimpleCortexRecord("AAAAACCAAATCATCTGTTGGAAATTTATTA", new int[] {0,1}, new String[] {"........", "a....C.."}));
        recs.add(new SimpleCortexRecord("AAACGATTTTCTATAGCTATGTAGTCATGCA", new int[] {1,1}, new String[] {"a.......", "a......."}));
        recs.add(new SimpleCortexRecord("AATTTATTCCAAATACTGCAAATACCCCAAA", new int[] {1,0}, new String[] {"a...A...", "........"}));
        recs.add(new SimpleCortexRecord("AATTTATTACAAATACTGCAAATACCCCAAA", new int[] {0,1}, new String[] {"........", "a...A..."}));
        recs.add(new SimpleCortexRecord("TACTGCAAATACCCCAAAACGATTTTCTATA", new int[] {1,1}, new String[] {"a.....G.", "a.....G."}));
        recs.add(new SimpleCortexRecord("TCTGTTGGAAATTTATTCCAAATACTGCAAA", new int[] {1,0}, new String[] {"a......T", "........"}));
        recs.add(new SimpleCortexRecord("CAAATCATCTGTTGGAAATTTATTACAAATA", new int[] {0,1}, new String[] {"........", ".c...C.."}));
        recs.add(new SimpleCortexRecord("ATCTGTTGGAAATTTATTACAAATACTGCAA", new int[] {0,1}, new String[] {"........", ".c..A..."}));
        recs.add(new SimpleCortexRecord("ATTTGCAGTATTTGTAATAAATTTCCAACAG", new int[] {0,1}, new String[] {"........", "...tA..."}));
        recs.add(new SimpleCortexRecord("ACCAAATCATCTGTTGGAAATTTATTCCAAA", new int[] {1,0}, new String[] {"a......T", "........"}));
        recs.add(new SimpleCortexRecord("CCCCAAAACGATTTTCTATAGCTATGTAGTC", new int[] {1,1}, new String[] {"a...A...", "a...A..."}));
        recs.add(new SimpleCortexRecord("AAATTTATTCCAAATACTGCAAATACCCCAA", new int[] {1,0}, new String[] {"..g.A...", "........"}));
        recs.add(new SimpleCortexRecord("ACTGCAAATACCCCAAAACGATTTTCTATAG", new int[] {1,1}, new String[] {"...t.C..", "...t.C.."}));
        recs.add(new SimpleCortexRecord("ATCGTTTTGGGGTATTTGCAGTATTTGTAAT", new int[] {0,1}, new String[] {"........", "a...A..."}));
        recs.add(new SimpleCortexRecord("AAAACCAAATCATCTGTTGGAAATTTATTAC", new int[] {0,1}, new String[] {"........", "a...A..."}));
        recs.add(new SimpleCortexRecord("GGGTATTTGCAGTATTTGGAATAAATTTCCA", new int[] {1,0}, new String[] {"..g.A...", "........"}));
        recs.add(new SimpleCortexRecord("AAATTTATTACAAATACTGCAAATACCCCAA", new int[] {0,1}, new String[] {"........", "..g.A..."}));
        recs.add(new SimpleCortexRecord("TATTTGCAGTATTTGTAATAAATTTCCAACA", new int[] {0,1}, new String[] {"........", "..g...G."}));
        recs.add(new SimpleCortexRecord("AAATCATCTGTTGGAAATTTATTCCAAATAC", new int[] {1,0}, new String[] {".c.....T", "........"}));
        recs.add(new SimpleCortexRecord("CATCTGTTGGAAATTTATTCCAAATACTGCA", new int[] {1,0}, new String[] {"...tA...", "........"}));
        recs.add(new SimpleCortexRecord("ATCATCTGTTGGAAATTTATTCCAAATACTG", new int[] {1,0}, new String[] {"a....C..", "........"}));
        recs.add(new SimpleCortexRecord("GGGTATTTGCAGTATTTGTAATAAATTTCCA", new int[] {0,1}, new String[] {"........", "..g.A..."}));
        recs.add(new SimpleCortexRecord("AAATACTGCAAATACCCCAAAACGATTTTCT", new int[] {1,1}, new String[] {".c..A...", ".c..A..."}));
        recs.add(new SimpleCortexRecord("ATACTGCAAATACCCCAAAACGATTTTCTAT", new int[] {1,1}, new String[] {"a...A...", "a...A..."}));
        recs.add(new SimpleCortexRecord("AAAACGATTTTCTATAGCTATGTAGTCATGC", new int[] {1,1}, new String[] {".c..A...", ".c..A..."}));
        recs.add(new SimpleCortexRecord("AATCGTTTTGGGGTATTTGCAGTATTTGTAA", new int[] {0,1}, new String[] {"........", "a......T"}));
        recs.add(new SimpleCortexRecord("CATCTGTTGGAAATTTATTACAAATACTGCA", new int[] {0,1}, new String[] {"........", "...tA..."}));
        recs.add(new SimpleCortexRecord("TATTCCAAATACTGCAAATACCCCAAAACGA", new int[] {1,0}, new String[] {"...t...T", "........"}));
        recs.add(new SimpleCortexRecord("AAAAACCAAATCATCTGTTGGAAATTTATTC", new int[] {1,0}, new String[] {"a....C..", "........"}));
        recs.add(new SimpleCortexRecord("ATAGCTATAGAAAATCGTTTTGGGGTATTTG", new int[] {1,1}, new String[] {".c...C..", ".c...C.."}));
        recs.add(new SimpleCortexRecord("GCAGTATTTGTAATAAATTTCCAACAGATGA", new int[] {0,1}, new String[] {"........", "...t...T"}));
        recs.add(new SimpleCortexRecord("CCCAAAACGATTTTCTATAGCTATGTAGTCA", new int[] {1,1}, new String[] {".c.....T", ".c.....T"}));
        recs.add(new SimpleCortexRecord("AATCGTTTTGGGGTATTTGCAGTATTTGGAA", new int[] {1,0}, new String[] {"a......T", "........"}));
        recs.add(new SimpleCortexRecord("AAAACCAAATCATCTGTTGGAAATTTATTCC", new int[] {1,0}, new String[] {"a...A...", "........"}));
        recs.add(new SimpleCortexRecord("AGCTATAGAAAATCGTTTTGGGGTATTTGCA", new int[] {1,1}, new String[] {"...t..G.", "...t..G."}));
        recs.add(new SimpleCortexRecord("CAAAACGATTTTCTATAGCTATGTAGTCATG", new int[] {1,1}, new String[] {".c...C..", ".c...C.."}));
        recs.add(new SimpleCortexRecord("TATTTGCAGTATTTGGAATAAATTTCCAACA", new int[] {1,0}, new String[] {"..g...G.", "........"}));
        recs.add(new SimpleCortexRecord("TATTACAAATACTGCAAATACCCCAAAACGA", new int[] {0,1}, new String[] {"........", "...t...T"}));
        recs.add(new SimpleCortexRecord("ATCTGTTGGAAATTTATTCCAAATACTGCAA", new int[] {1,0}, new String[] {".c..A...", "........"}));
        recs.add(new SimpleCortexRecord("AATACTGCAAATACCCCAAAACGATTTTCTA", new int[] {1,1}, new String[] {"a......T", "a......T"}));
        recs.add(new SimpleCortexRecord("ACCCCAAAACGATTTTCTATAGCTATGTAGT", new int[] {1,1}, new String[] {"...t.C..", "...t.C.."}));
        recs.add(new SimpleCortexRecord("CTACATAGCTATAGAAAATCGTTTTGGGGTA", new int[] {1,1}, new String[] {"a......T", "a......T"}));
        recs.add(new SimpleCortexRecord("ATTTGGAATAAATTTCCAACAGATGATTTGG", new int[] {1,0}, new String[] {"...t...T", "........"}));
        recs.add(new SimpleCortexRecord("ACCAAATCATCTGTTGGAAATTTATTACAAA", new int[] {0,1}, new String[] {"........", "a......T"}));
        recs.add(new SimpleCortexRecord("AATCATCTGTTGGAAATTTATTCCAAATACT", new int[] {1,0}, new String[] {"a.....G.", "........"}));
        recs.add(new SimpleCortexRecord("ATGACTACATAGCTATAGAAAATCGTTTTGG", new int[] {1,1}, new String[] {".c....G.", ".c....G."}));
        recs.add(new SimpleCortexRecord("AATCATCTGTTGGAAATTTATTACAAATACT", new int[] {0,1}, new String[] {"........", "a.....G."}));
        recs.add(new SimpleCortexRecord("CGTTTTGGGGTATTTGCAGTATTTGTAATAA", new int[] {0,1}, new String[] {"........", "...tA..."}));
        recs.add(new SimpleCortexRecord("GGAAATTTATTACAAATACTGCAAATACCCC", new int[] {0,1}, new String[] {"........", "...tA..."}));
        recs.add(new SimpleCortexRecord("GAAATTTATTACAAATACTGCAAATACCCCA", new int[] {0,1}, new String[] {"........", "..g.A..."}));
        recs.add(new SimpleCortexRecord("GGTATTTGCAGTATTTGGAATAAATTTCCAA", new int[] {1,0}, new String[] {"..g..C..", "........"}));
        recs.add(new SimpleCortexRecord("GTTTTGGGGTATTTGCAGTATTTGTAATAAA", new int[] {0,1}, new String[] {"........", ".c.....T"}));
        recs.add(new SimpleCortexRecord("AACCAAATCATCTGTTGGAAATTTATTCCAA", new int[] {1,0}, new String[] {"a...A...", "........"}));

        CortexGraph cg = new CortexGraph("testdata/smallgraph.ctx");

        int index = 0;
        for (CortexRecord cr : cg) {
            SimpleCortexRecord scr = recs.get(index);

            Assert.assertEquals(scr.equals(cr), true, "Cortex record says '" + cr + "' but test record says '" + scr + "'");

            index++;
        }
    }
}
