package uk.ac.ox.well.indiana.utils.assembly;

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.util.*;

public class CortexGraphWalkerTest {
    private CortexGraph cg;
    private Map<String, CortexRecord> records;
    private CortexGraphWalker cgw;

    @BeforeClass
    public void setup() {
        records = new HashMap<String, CortexRecord>();

        cg = new CortexGraph("testdata/test_gene_for_sn_reconstruction.ctx");
        for (CortexRecord cr : cg) {
            String kmer = cr.getKmerString();
            records.put(kmer, cr);
        }

        cgw = new CortexGraphWalker(records);
    }

    @Test
    public void testSupernodeReconstruction() {
        String referenceKmer = SequenceUtils.getAlphanumericallyLowestOrientation("GTTTCCACAGCGGCCGTTTCCACGAAGGCGG");

        String supernode = cgw.getSupernode(0, referenceKmer);

        Assert.assertEquals("CCGCCTTCGTGGAAACGGCCGCTGTGGAAACTTTTTTCTTATGGCATAAGTATAAACAAGAGAAGAAGAAACCAAAAAATGAAGTGGGAGGTGCAGCGGGAGTACTACAAA", supernode);
    }

    @Test
    public void testReferenceGuidedSupernodeReconstruction() {
        Set<String> referenceKmers = new HashSet<String>();
        referenceKmers.add(SequenceUtils.getAlphanumericallyLowestOrientation("GTTTCCACAGCGGCCGTTTCCACGAAGGCGG"));
        referenceKmers.add(SequenceUtils.getAlphanumericallyLowestOrientation("GTTTGTAGTACTCCCGCTGCACCTCCCACTT"));
        referenceKmers.add(SequenceUtils.getAlphanumericallyLowestOrientation("AAGGCGGTAAGGAGGAGCGCGTCACCTTTTG"));

        Collection<String> supernodes = cgw.getReferenceGuidedSupernodes(0, referenceKmers);

        Assert.assertEquals(1, supernodes.size());
        Assert.assertEquals("CCGCCGATGGTTTGTAGTACTCCCGCTGCACCTCCCACTTCATTTTTTGGTTTCTTCTTCTCTTGTTTATACTTATGCCATAAGAAAAAAGTTTCCACAGCGGCCGTTTCCACGAAGGCGGTAAGGAGGAGCGCGTCACCTTTTGGTGACTGCGACGGAGTTCCCGCCTGAGTGGCTTCA", supernodes.iterator().next());
    }

    @Test
    public void testMultipleReferenceGuidedSupernodeReconstruction() {
        Set<String> referenceKmers = new HashSet<String>();
        referenceKmers.add(SequenceUtils.getAlphanumericallyLowestOrientation("GTTTCCACAGCGGCCGTTTCCACGAAGGCGG"));
        referenceKmers.add(SequenceUtils.getAlphanumericallyLowestOrientation("GTTTGTAGTACTCCCGCTGCACCTCCCACTT"));
        referenceKmers.add(SequenceUtils.getAlphanumericallyLowestOrientation("AAGGCGGTAAGGAGGAGCGCGTCACCTTTTG"));
        referenceKmers.add(SequenceUtils.getAlphanumericallyLowestOrientation("TGCTGTTTCTAGGCGGCCATTTGTAGTACTC"));

        Collection<String> supernodes = cgw.getReferenceGuidedSupernodes(0, referenceKmers);

        Assert.assertEquals(2, supernodes.size());
        Assert.assertEquals(true, supernodes.contains("GAAGTGGGAGGTGCAGCGGGAGTACTACAAATGGCCGCCTAGAAACAGCA"));
        Assert.assertEquals(true, supernodes.contains("CCGCCGATGGTTTGTAGTACTCCCGCTGCACCTCCCACTTCATTTTTTGGTTTCTTCTTCTCTTGTTTATACTTATGCCATAAGAAAAAAGTTTCCACAGCGGCCGTTTCCACGAAGGCGGTAAGGAGGAGCGCGTCACCTTTTGGTGACTGCGACGGAGTTCCCGCCTGAGTGGCTTCA"));
    }
}
