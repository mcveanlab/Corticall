package uk.ac.ox.well.indiana.utils.sequence;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class SequenceUtilsTest {
    @Test
    public void complementTest() {
        byte[] trials = new byte[] { 'A', 'C', 'G', 'T', 'N', '.', 'a', 'c', 'g', 't' };
        byte[] exp    = new byte[] { 'T', 'G', 'C', 'A', 'N', '.', 't', 'g', 'c', 'a' };

        for (int i = 0; i < trials.length; i++) {
            Assert.assertEquals(SequenceUtils.complement(trials[i]), exp[i]);
        }
    }

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
        for (int i = 0; i < 10000; i++) {
            for (Integer k : Arrays.asList(21, 31, 41, 51)) {
                String fw = new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(k));
                String rc = SequenceUtils.reverseComplement(fw);

                byte[] expected = fw.compareTo(rc) < 0 ? fw.getBytes() : rc.getBytes();
                byte[] computed = SequenceUtils.alphanumericallyLowestOrientation(fw.getBytes());

                Assert.assertEquals(expected, computed);
            }
        }
    }

    @Test
    public void computeN50() {
        List<String> sequences = new ArrayList<>();
        sequences.add("AAGCTTA");
        sequences.add("TTGA");
        sequences.add("AAC");
        sequences.add("TT");
        sequences.add("AA");
        sequences.add("C");
        sequences.add("G");

        int n50 = SequenceUtils.computeN50Value(sequences);

        Assert.assertEquals(4, n50);

        List<CortexKmer> csequences = new ArrayList<>();
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

            Set<String> cdsIds = new HashSet<>();
            Set<String> trIds = new HashSet<>();

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

    @Test
    public void testRegex() {
        String template = "...99999999999999999999999999999999999999999999999999_____________";

        String pattern = "\\.+_*([A-Za-z0-9])\\1+.+";
        Pattern p = Pattern.compile(pattern);
        Matcher m = p.matcher(template);

        System.out.println(template);
        System.out.println(m);
        System.out.println(m.matches());
        System.out.println(m.groupCount());
        for (int i = 0; i <= m.groupCount(); i++) {
            System.out.println("  group=" + m.group(i));
        }
    }
}
