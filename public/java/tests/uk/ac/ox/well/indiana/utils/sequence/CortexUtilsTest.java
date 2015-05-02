package uk.ac.ox.well.indiana.utils.sequence;

import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexJunctionsRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksMap;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksRecord;
import java.util.List;

public class CortexUtilsTest {
    @Test
    public void testKmersInLink() {
        CortexGraph cg = new CortexGraph("testdata/PG0051-C.ERR019061.chr1.infer.sorted.ctx");
        CortexLinksMap ctp = new CortexLinksMap("testdata/PG0051-C.ERR019061.chr1.se.ctp");

        for (CortexKmer ock : ctp.keySet()) {
            for (int rc = 0; rc <= 1; rc++) {
                String sk = ock.getKmerAsString();
                if (rc == 1) { sk = SequenceUtils.reverseComplement(sk); }

                CortexLinksRecord clr = ctp.get(sk);

                for (CortexJunctionsRecord cjr : clr.getJunctions()) {
                    List<String> kils = CortexUtils.getKmersInLinkByNavigation(cg, sk, cjr);
                    List<String> expectedKmers = CortexUtils.getKmersInLinkFromSeq(cg, sk, cjr);

                    for (int i = 0; i < kils.size(); i++) {
                        String expectedKmer = expectedKmers.get(i);
                        String seenKmer = kils.get(i);

                        Assert.assertEquals(seenKmer, expectedKmer);
                    }
                }
            }
        }
    }

    @Test(enabled = false)
    public void decodeBinaryKmer() {
        CortexGraph cg = new CortexGraph("testdata/PG0051-C.ERR019061.chr1.infer.sorted.ctx");

        for (CortexRecord cr : cg) {
            System.out.println(cr);
        }
    }
}
