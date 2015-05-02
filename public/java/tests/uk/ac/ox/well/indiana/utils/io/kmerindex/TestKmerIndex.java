package uk.ac.ox.well.indiana.utils.io.kmerindex;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Interval;
import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.fasta.MultiIndexedFastaSequenceFile;

import java.io.File;
import java.io.FileNotFoundException;

public class TestKmerIndex {
    private final File testChr = new File("testdata/Pf3D7_01_v3.fasta");

    @Test
    public void testWriteIndex() {
        int kmerSize = 47;
        int numKmersWritten = KmerIndex.writeIndex(testChr, kmerSize);
        int numKmersFound = 0;

        MultiIndexedFastaSequenceFile iref;
        try {
            iref = new MultiIndexedFastaSequenceFile(testChr);
        } catch (FileNotFoundException e) {
            throw new IndianaException("File not found: " + e);
        }

        ReferenceSequenceFile ref = new FastaSequenceFile(testChr, true);
        ReferenceSequence rseq;
        while ((rseq = ref.nextSequence()) != null) {
            String seq = new String(rseq.getBases());

            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                String kmer = seq.substring(i, i + kmerSize);

                Interval interval = iref.find(kmer);
                if (interval != null) {
                    Assert.assertEquals(interval.getStart(), i);
                    numKmersFound++;
                }
            }
        }

        Assert.assertEquals(numKmersFound, numKmersWritten);
    }
}
