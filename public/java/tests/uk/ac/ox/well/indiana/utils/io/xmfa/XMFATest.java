package uk.ac.ox.well.indiana.utils.io.xmfa;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

public class XMFATest {
    @Test
    public void loadXMFAFileTest() {
        FastaSequenceFile testRecordFile = new FastaSequenceFile(new File("testdata/pf.subset.alignment.fasta"), true);

        Map<String, String> testRecords = new HashMap<String, String>();
        ReferenceSequence rseq;
        while ((rseq = testRecordFile.nextSequence()) != null) {
            String name = rseq.getName();
            String seq = new String(rseq.getBases());

            testRecords.put(name, seq);
        }

        XMFASequenceFile xf = new XMFASequenceFile("testdata/pf.alignment");

        Assert.assertEquals(21856, xf.getNumRecords());

        XMFARecord xr = xf.iterator().next();
        for (String name : xr.keySet()) {
            String seq = new String(xr.get(name).getBases());

            Assert.assertEquals(testRecords.get(name), seq);
        }
    }
}
