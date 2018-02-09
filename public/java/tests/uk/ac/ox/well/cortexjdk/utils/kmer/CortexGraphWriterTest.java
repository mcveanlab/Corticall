package uk.ac.ox.well.cortexjdk.utils.kmer;

import com.google.common.collect.Lists;
import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.cortexjdk.utils.assembler.TempGraphAssembler;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexColor;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraphWriter;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class CortexGraphWriterTest {
    @Test
    public void writeSmallGraphTest() {
        try {
            Map<String, Collection<String>> haplotypes = new HashMap<>();
            haplotypes.put("mom", Collections.singletonList(new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(100))));
            haplotypes.put("dad", Collections.singletonList(new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(100))));
            haplotypes.put("kid", Collections.singletonList(new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(100))));
            CortexGraph cg1 = TempGraphAssembler.buildGraph(haplotypes, 5);

            File tempFile = File.createTempFile("smallgraph-copy", ".ctx");
            tempFile.deleteOnExit();

            CortexGraphWriter cgw = new CortexGraphWriter(tempFile);
            cgw.setHeader(cg1.getHeader());

            for (CortexRecord cr : cg1) {
                cgw.addRecord(cr);
            }

            cgw.close();

            CortexGraph cg2 = new CortexGraph(tempFile);

            Assert.assertEquals(cg1.getVersion(), cg2.getVersion());
            Assert.assertEquals(cg1.getKmerBits(), cg2.getKmerBits());
            Assert.assertEquals(cg1.getKmerSize(), cg2.getKmerSize());
            Assert.assertEquals(cg1.getNumColors(), cg2.getNumColors());
            Assert.assertEquals(cg1.getNumRecords(), cg2.getNumRecords());

            for (int c = 0; c < cg1.getNumColors(); c++) {
                CortexColor cc1 = cg1.getColor(c);
                CortexColor cc2 = cg2.getColor(c);

                Assert.assertEquals(cc1.getCleanedAgainstGraphName(), cc2.getCleanedAgainstGraphName());
                Assert.assertEquals(cc1.getLowCovKmerThreshold(), cc2.getLowCovKmerThreshold());
                Assert.assertEquals(cc1.getLowCovSupernodesThreshold(), cc2.getLowCovSupernodesThreshold());
                Assert.assertEquals(cc1.getMeanReadLength(), cc2.getMeanReadLength());
                Assert.assertEquals(cc1.getTotalSequence(), cc2.getTotalSequence());
                Assert.assertEquals(cc1.getSampleName(), cc2.getSampleName());
                Assert.assertEquals(cc1.isCleanedAgainstGraph(), cc2.isCleanedAgainstGraph());
                Assert.assertEquals(cc1.isLowCovgKmersRemoved(), cc2.isLowCovgKmersRemoved());
                Assert.assertEquals(cc1.isLowCovgSupernodesRemoved(), cc2.isLowCovgSupernodesRemoved());
                Assert.assertEquals(cc1.isTipClippingApplied(), cc2.isTipClippingApplied());
            }

            Iterator<CortexRecord> it1 = cg1.iterator();
            Iterator<CortexRecord> it2 = cg2.iterator();

            while (it1.hasNext() && it2.hasNext()) {
                CortexRecord cr1 = it1.next();
                CortexRecord cr2 = it2.next();

                Assert.assertEquals(cr1, cr2);
            }

            Assert.assertEquals(it1.hasNext(), it2.hasNext());
        } catch (IOException e) {
            throw new CortexJDKException("Could not create temp file", e);
        }
    }
}
