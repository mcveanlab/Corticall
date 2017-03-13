package uk.ac.ox.well.indiana.utils.io.cortex.graph;

import com.google.common.base.Joiner;
import org.apache.commons.math3.util.Pair;
import org.junit.Test;
import org.testng.Assert;

import java.util.*;

/**
 * Created by kiran on 13/03/2017.
 */
public class CortexIndexTest {
    @Test
    public void validateIndex() {
        CortexGraph graph = new CortexGraph("testdata/laverania.ctx", false);
        CortexIndex index = new CortexIndex("testdata/laverania.ctx.idx", graph);

        Map<String, Pair<Long, Long>> mappedIndex = new HashMap<>();

        for (Pair<String, Pair<Long, Long>> p : index.getIndex()) {
            mappedIndex.put(p.getFirst(), p.getSecond());
        }

        for (long i = 0; i < graph.getNumRecords(); i++) {
            CortexRecord cr = graph.getRecord(i);
            String sk = cr.getKmerAsString();

            if (mappedIndex.containsKey(sk)) {
                Assert.assertEquals(mappedIndex.get(sk).getFirst().longValue(), i);
            }
        }
    }

    @Test
    public void retreiveWithIndex() {
        CortexGraph gIgnoreIndex = new CortexGraph("testdata/laverania.ctx", true);
        CortexGraph gUseIndex = new CortexGraph("testdata/laverania.ctx", false);

        List<CortexKmer> kmers = new ArrayList<>();
        for (CortexRecord cr : gIgnoreIndex) {
            kmers.add(cr.getCortexKmer());
        }

        kmers.sort(new Comparator<CortexKmer>() {
            @Override
            public int compare(CortexKmer o1, CortexKmer o2) {
                return Integer.compare(o1.hashCode(), o2.hashCode());
            }
        });

        CortexRecord cIgnoreIndex = gIgnoreIndex.findRecord(kmers.get(0));
        CortexRecord cUseIndex = gUseIndex.findRecord(kmers.get(0));
        Assert.assertEquals(cIgnoreIndex, cUseIndex);

        cIgnoreIndex = gIgnoreIndex.findRecord(kmers.get(kmers.size() - 1));
        cUseIndex = gUseIndex.findRecord(kmers.get(kmers.size() - 1));
        Assert.assertEquals(cIgnoreIndex, cUseIndex);

        for (int i = 0; i < kmers.size(); i += 100) {
            cIgnoreIndex = gIgnoreIndex.findRecord(kmers.get(i));
            cUseIndex = gUseIndex.findRecord(kmers.get(i));

            Assert.assertEquals(cIgnoreIndex, cUseIndex);
        }
    }

    @Test
    public void retreiveWithIndexSmallFile() {
        CortexGraph gIgnoreIndex = new CortexGraph("testdata/laverania_contig_1.ctx", true);
        CortexGraph gUseIndex = new CortexGraph("testdata/laverania_contig_1.ctx", false);

        List<CortexKmer> kmers = new ArrayList<>();
        for (CortexRecord cr : gIgnoreIndex) {
            kmers.add(cr.getCortexKmer());
        }

        kmers.sort(new Comparator<CortexKmer>() {
            @Override
            public int compare(CortexKmer o1, CortexKmer o2) {
                return Integer.compare(o1.hashCode(), o2.hashCode());
            }
        });

        CortexRecord cIgnoreIndex = gIgnoreIndex.findRecord(kmers.get(0));
        CortexRecord cUseIndex = gUseIndex.findRecord(kmers.get(0));
        Assert.assertEquals(cIgnoreIndex, cUseIndex);

        cIgnoreIndex = gIgnoreIndex.findRecord(kmers.get(kmers.size() - 1));
        cUseIndex = gUseIndex.findRecord(kmers.get(kmers.size() - 1));
        Assert.assertEquals(cIgnoreIndex, cUseIndex);

        for (int i = 0; i < kmers.size(); i++) {
            cIgnoreIndex = gIgnoreIndex.findRecord(kmers.get(i));
            cUseIndex = gUseIndex.findRecord(kmers.get(i));

            Assert.assertEquals(cIgnoreIndex, cUseIndex);
        }
    }
}
