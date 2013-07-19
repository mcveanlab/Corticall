package uk.ac.ox.well.indiana.utils.assembly;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import org.jgrapht.DirectedGraph;
import org.jgrapht.ext.DOTExporter;
import org.jgrapht.graph.DefaultEdge;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexMap;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class CortexGraphWalkerTest {
    private CortexGraph cg;
    private Map<String, CortexRecord> records;
    private CortexGraphWalker cgw;

    private String cg14Seq;
    private CortexMap cg14Map;
    private CortexGraphWalker cg14Walker;

    private DirectedGraph<CortexKmer, DefaultEdge> gfw;
    private DirectedGraph<CortexKmer, DefaultEdge> grc;

    private CortexKmer ckfw;
    private CortexKmer ckrc;

    @BeforeClass
    public void setup() {
        records = new HashMap<String, CortexRecord>();

        cg = new CortexGraph("testdata/test_gene_for_sn_reconstruction.ctx");
        for (CortexRecord cr : cg) {
            String kmer = cr.getKmerAsString();
            records.put(kmer, cr);
        }

        cgw = new CortexGraphWalker(new CortexMap(cg));

        FastaSequenceFile cg14Fasta = new FastaSequenceFile(new File("testdata/Pf3D7_14_v3.fasta"), true);
        cg14Seq = new String(cg14Fasta.nextSequence().getBases());
        cg14Map = new CortexMap(new File("testdata/Pf3D7_14_v3.ctx"));
        cg14Walker = new CortexGraphWalker(cg14Map);

        String fw = "AACGGCCGCTGTGGAAACTTTTTTCTTATGG";
        String rc = SequenceUtils.reverseComplement(fw);

        ckfw = new CortexKmer(fw);
        ckrc = new CortexKmer(rc);
    }

    @Test
    public void testContigReconstruction() {
        CortexKmer panelKmer = new CortexKmer("GTTTCCACAGCGGCCGTTTCCACGAAGGCGG");

        CortexKmer contig = cgw.buildContig(0, panelKmer);

        Assert.assertEquals("CCGCCTTCGTGGAAACGGCCGCTGTGGAAACTTTTTTCTTATGGCATAAGTATAAACAAGAGAAGAAGAAACCAAAAAATGAAGTGGGAGGTGCAGCGGGAGTACTACAAA", contig);
    }

    @Test
    public void testRandomReconstructedContigsMatchOriginExactly() {
        Random random = new Random();

        Set<CortexKmer> panel = new HashSet<CortexKmer>();

        for (CortexRecord cr : cg14Map.values()) {
            if (cr.getCoverage(0) == 1) {
                if (random.nextBoolean()) {
                    panel.add(cr.getKmer());
                }
            }

            if (panel.size() >= 100) {
                break;
            }
        }

        Map<CortexKmer, Set<CortexKmer>> contigs = cg14Walker.buildContigs(0, panel);

        for (CortexKmer contig : contigs.keySet()) {
            String fw = contig.getKmerAsString();
            String rc = SequenceUtils.reverseComplement(fw);

            Assert.assertTrue(cg14Seq.contains(fw) || cg14Seq.contains(rc));
        }
    }

    @Test
    public void testLocalGraphContainsVerticesAndEdgesFromFastaFile() {
        gfw = cgw.buildLocalGraph(0, ckfw, 1);
        grc = cgw.buildLocalGraph(0, ckrc, 1);

        int kmerSize = gfw.vertexSet().iterator().next().length();

        FastaSequenceFile fa = new FastaSequenceFile(new File("testdata/test_gene_for_sn_reconstruction.fasta"), true);
        ReferenceSequence seq;
        while ((seq = fa.nextSequence()) != null) {
            String bases = new String(seq.getBases());

            for (int i = 0; i < bases.length() - kmerSize; i++) {
                String kmer0 = bases.substring(i, i + kmerSize);
                CortexKmer ck0 = new CortexKmer(kmer0);

                String kmer1 = bases.substring(i + 1, i + 1 + kmerSize);
                CortexKmer ck1 = new CortexKmer(kmer1);

                Assert.assertTrue(gfw.containsVertex(ck0), "Fw graph does not contain vertex '" + ck0 + "'");
                Assert.assertTrue(gfw.containsVertex(ck1), "Fw graph does not contain vertex '" + ck1 + "'");
                Assert.assertTrue(gfw.containsEdge(ck0, ck1) || gfw.containsEdge(ck1, ck0), "Fw graph does not contain edge ['" + ck0 + "', '" + ck1 + "']");

                Assert.assertTrue(grc.containsVertex(ck0), "Rc graph does not contain vertex '" + ck0 + "'");
                Assert.assertTrue(grc.containsVertex(ck1), "Rc graph does not contain vertex '" + ck1 + "'");
                Assert.assertTrue(grc.containsEdge(ck0, ck1) || grc.containsEdge(ck1, ck0), "Rc graph does not contain edge ['" + ck0 + "', '" + ck1 + "']");
            }
        }

        for (CortexKmer ck : gfw.vertexSet()) {
            Assert.assertTrue(grc.containsVertex(ck), "Fw graph contains a vertex that the rc graph does not ('" + ck + "')");
        }
    }

    @Test
    public void testLocalGraphContainsExpectedForks() {
        CortexKmer aPanelKmer = new CortexKmer("CCGCCTTCGTGGAAACGGCCGCTGTGGAAAC");
        CortexKmer aIncomingKmer1 = new CortexKmer("CGCCTTCGTGGAAACGGCCGCTGTGGAAACT");
        CortexKmer aOutgoingKmer1 = new CortexKmer("ACCGCCTTCGTGGAAACGGCCGCTGTGGAAA");
        CortexKmer aOutgoingKmer2 = new CortexKmer("GCCGCCTTCGTGGAAACGGCCGCTGTGGAAA");

        Assert.assertTrue(gfw.containsVertex(aPanelKmer));
        Assert.assertTrue(gfw.containsEdge(aIncomingKmer1, aPanelKmer));
        Assert.assertTrue(gfw.containsEdge(aPanelKmer, aOutgoingKmer1));
        Assert.assertTrue(gfw.containsEdge(aPanelKmer, aOutgoingKmer2));
        Assert.assertEquals(gfw.inDegreeOf(aPanelKmer), 1);
        Assert.assertEquals(gfw.outDegreeOf(aPanelKmer), 2);

        CortexKmer bPanelKmer = new CortexKmer("GAAGTGGGAGGTGCAGCGGGAGTACTACAAA");
        CortexKmer bIncomingKmer1 = new CortexKmer("TGAAGTGGGAGGTGCAGCGGGAGTACTACAA");
        CortexKmer bOutgoingKmer1 = new CortexKmer("AAGTGGGAGGTGCAGCGGGAGTACTACAAAC");
        CortexKmer bOutgoingKmer2 = new CortexKmer("AAGTGGGAGGTGCAGCGGGAGTACTACAAAT");

        Assert.assertTrue(gfw.containsVertex(bPanelKmer));
        Assert.assertTrue(gfw.containsEdge(bIncomingKmer1, bPanelKmer));
        Assert.assertTrue(gfw.containsEdge(bPanelKmer, bOutgoingKmer1));
        Assert.assertTrue(gfw.containsEdge(bPanelKmer, bOutgoingKmer2));
        Assert.assertEquals(gfw.inDegreeOf(bPanelKmer), 1);
        Assert.assertEquals(gfw.outDegreeOf(bPanelKmer), 2);
    }

    @Test
    public void testMultiSampleGraph() {
        System.out.println("testTwoSampleGraph()");

        CortexMap cm = new CortexMap("testdata/samples_1_and_2.ctx");
        CortexKmer ck = new CortexKmer("GGATAAATCAAGTATTGCTAACAAAATTGAA");
        CortexGraphWalker gw = new CortexGraphWalker(cm);

        DirectedGraph<CortexKmer, DefaultEdge> g = gw.buildLocalGraph(ck, 2);

        DOTExporter<CortexKmer, DefaultEdge> exporter = new DOTExporter<CortexKmer, DefaultEdge>(new CortexKmerIDProvider(), new CortexKmerLabelProvider(), null, new CortexKmerAttributeProvider(ck), null);
        try {
            exporter.export(new FileWriter("test10.dot"), g);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
