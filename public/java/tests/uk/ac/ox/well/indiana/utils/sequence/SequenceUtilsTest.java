package uk.ac.ox.well.indiana.utils.sequence;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.*;
import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class SequenceUtilsTest {
    @Test
    public void reverseComplementTest() {
        byte[] sequence   = "TACTGACTTTTCTCGCTATTCGTATGCATG".getBytes();
        byte[] expectedRC = "CATGCATACGAATAGCGAGAAAAGTCAGTA".getBytes();

        byte[] computedRC = SequenceUtils.reverseComplement(sequence);

        Assert.assertEquals(expectedRC, computedRC);
    }

    @Test
    public void reverseComplementTestWithNs() {
        byte[] sequence   = "NACTGACTTTTCTCGCTATTCGTATGCATG".getBytes();
        byte[] expectedRC = "CATGCATACGAATAGCGAGAAAAGTCAGTN".getBytes();

        byte[] computedRC = SequenceUtils.reverseComplement(sequence);

        Assert.assertEquals(expectedRC, computedRC);
    }

    @Test
    public void reverseComplementWithLowercaseBases() {
        byte[] sequence   = "NACTGACTTTTCTCGCTATTCGTATGCATg".getBytes();
        byte[] expectedRC = "cATGCATACGAATAGCGAGAAAAGTCAGTN".getBytes();

        byte[] computedRC = SequenceUtils.reverseComplement(sequence);

        Assert.assertEquals(expectedRC, computedRC);
    }

    @Test
    public void alphanumericallyLowestOrientation() {
        byte[] sequence            = "TACTGACTTTTCTCGCTATTCGTATGCATG".getBytes();
        byte[] expectedOrientation = "CATGCATACGAATAGCGAGAAAAGTCAGTA".getBytes();

        byte[] computedOrientation = SequenceUtils.alphanumericallyLowestOrientation(sequence);

        Assert.assertEquals(expectedOrientation, computedOrientation);
    }

    @Test
    public void computeN50() {
        List<String> sequences = new ArrayList<String>();
        sequences.add("AAGCTTA");
        sequences.add("TTGA");
        sequences.add("AAC");
        sequences.add("TT");
        sequences.add("AA");
        sequences.add("C");
        sequences.add("G");

        int n50 = SequenceUtils.computeN50Value(sequences);

        Assert.assertEquals(4, n50);

        List<CortexKmer> csequences = new ArrayList<CortexKmer>();
        csequences.add(new CortexKmer("AAGCTTA"));
        csequences.add(new CortexKmer("TTGA"));
        csequences.add(new CortexKmer("AAC"));
        csequences.add(new CortexKmer("TT"));
        csequences.add(new CortexKmer("AA"));
        csequences.add(new CortexKmer("C"));
        csequences.add(new CortexKmer("G"));

        int cn50 = SequenceUtils.computeN50Value(csequences);

        Assert.assertEquals(4, cn50);
    }

    @Test(enabled = false)
    public void testTranslateCodingSequence() {
        String cds1 = "ATGGTGACTGGTAGTGGTGGTGAGGATAAGTATAAAAGTGCCAAAAATGCCAAGGAACTT";
        String tr1  = "MVTGSGGEDKYKSAKNAKEL";

        Assert.assertEquals(tr1, SequenceUtils.translateCodingSequence(cds1));
    }

    @Test
    public void testExtractCodingSequenceAndTranslation() {
        File fastaFile = new File("testdata/Pf3D7_14_v3.fasta");
        File cdsFile = new File("testdata/Pf3D7_14_v3.cds.fasta");
        File trFile = new File("testdata/Pf3D7_14_v3.tr.fasta");
        File gffFile = new File("testdata/Pf3D7_14_v3.gff");

        try {
            IndexedFastaSequenceFile fastaRef = new IndexedFastaSequenceFile(fastaFile);
            IndexedFastaSequenceFile cdsRef = new IndexedFastaSequenceFile(cdsFile);
            IndexedFastaSequenceFile trRef = new IndexedFastaSequenceFile(trFile);

            Set<String> cdsIds = new HashSet<String>();
            Set<String> trIds = new HashSet<String>();

            ReferenceSequence r;
            while ((r = cdsRef.nextSequence()) != null) {
                String[] names = r.getName().split("\\s+");
                String name = names[0];

                cdsIds.add(name);
            }

            while ((r = trRef.nextSequence()) != null) {
                String[] names = r.getName().split("\\s+");
                String name = names[0];

                trIds.add(name);
            }

            GFF3 gff = new GFF3(gffFile);

            for (GFF3Record gr : gff) {
                if (gr.getType().equals("gene")) {
                    String id = gr.getAttribute("ID");

                    if (trIds.contains(id) && cdsIds.contains(id)) {
                        String cdsExpected = new String(cdsRef.getSequence(id).getBases());
                        String trExpected = new String(trRef.getSequence(id).getBases());

                        Collection<GFF3Record> exons = GFF3.getType("exon", gff.getChildren(gr));

                        String cds = SequenceUtils.extractCodingSequence(exons, fastaRef);
                        String tr = SequenceUtils.translateCodingSequence(cds);

                        Assert.assertEquals(cds, cdsExpected);
                        Assert.assertEquals(tr, trExpected);
                    }
                }
            }
        } catch (FileNotFoundException e) {
            throw new IndianaException("Error in opening test file '" + fastaFile.getAbsolutePath() + "'", e);
        }
    }
}
