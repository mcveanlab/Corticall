package uk.ac.ox.well.cortexjdk.utils.alignment.pairwise;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by kiran on 29/09/2017.
 */
public class TestBwaAligner {
    @Test
    public void testAlignment() {
        String ref = "tests/two_short_contigs.fa";
        Map<String, String> expectedAlignments = new HashMap<>();
        expectedAlignments.put("TTCTGTATCGTATGCTCTGAATAAAAATCGTGGCCCTATTTCGTATAGT", "unknown\t0\t1\t0\t60\t49M\t*\t0\t0\tTTCTGTATCGTATGCTCTGAATAAAAATCGTGGCCCTATTTCGTATAGT\t*\tNM:i:0");
        expectedAlignments.put("GGGCCGCGCCTATTATGGGCTTCTCTCTGAGTACTGGTCATGTAGTTGCTGTAGTCGTAGTGTCGTGGCCCCCCAGT", "unknown\t0\t2\t0\t60\t77M\t*\t0\t0\tGGGCCGCGCCTATTATGGGCTTCTCTCTGAGTACTGGTCATGTAGTTGCTGTAGTCGTAGTGTCGTGGCCCCCCAGT\t*\tNM:i:0");

        BwaAligner ba = new BwaAligner(ref);

        FastaSequenceFile fa = new FastaSequenceFile(new File(ref), true);
        ReferenceSequence rseq;
        while ((rseq = fa.nextSequence()) != null) {
            List<SAMRecord> alignments = ba.align(rseq.getBaseString());

            for (SAMRecord sr : alignments) {
                Assert.assertEquals(expectedAlignments.get(sr.getReadString()), sr.getSAMString().trim());
            }
        }

        ba.close();
    }
}
