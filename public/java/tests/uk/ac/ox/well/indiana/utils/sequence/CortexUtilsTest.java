package uk.ac.ox.well.indiana.utils.sequence;

import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexJunctionsRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksMap;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksRecord;

import java.util.*;

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

    @Test
    public void testNextKmers() {
        CortexGraph clean = new CortexGraph("testdata/smallgraph.sorted.ctx");

        Set<String> kmers = new HashSet<>();
        for (CortexRecord cr : clean) {
            kmers.add(cr.getKmerAsString());
        }

        for (String kmer : kmers) {
            Map<Integer, Set<String>> nkm = new HashMap<>();
            nkm.put(0, CortexUtils.getNextKmers(clean, null, kmer, 0));
            nkm.put(1, CortexUtils.getNextKmers(clean, null, kmer, 1));

            Map<Integer, Set<String>> nka = CortexUtils.getNextKmers(clean, null, kmer);

            Assert.assertEquals(nkm, nka);
        }
    }

    @Test
    public void textPrevKmers() {
        CortexGraph clean = new CortexGraph("testdata/smallgraph.sorted.ctx");

        Set<String> kmers = new HashSet<>();
        for (CortexRecord cr : clean) {
            kmers.add(cr.getKmerAsString());
        }

        for (String kmer : kmers) {
            Map<Integer, Set<String>> pkm = new HashMap<>();
            pkm.put(0, CortexUtils.getPrevKmers(clean, null, kmer, 0));
            pkm.put(1, CortexUtils.getPrevKmers(clean, null, kmer, 1));

            Map<Integer, Set<String>> pka = CortexUtils.getPrevKmers(clean, null, kmer);

            Assert.assertEquals(pkm, pka);
        }
    }
}
