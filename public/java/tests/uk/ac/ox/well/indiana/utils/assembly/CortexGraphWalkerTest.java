package uk.ac.ox.well.indiana.utils.assembly;

import net.sf.picard.reference.FastaSequenceFile;
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
    public void testBuildLocalGraph() {
        CortexKmer panelKmer = new CortexKmer("TCTTCTCTTGTTTATACTTATGCCATAAGAA");

        DirectedGraph<CortexKmer, DefaultEdge> graph = cgw.buildLocalGraph(0, panelKmer, 1);

        DOTExporter<CortexKmer, DefaultEdge> exporter = new DOTExporter<CortexKmer, DefaultEdge>(new CortexKmerNameProvider(), new CortexKmerNameProvider(), null);

        try {
            exporter.export(new FileWriter("initial-graph.dot"), graph);
        } catch (IOException e) {
            e.printStackTrace();
        }

        //System.out.println(graph);
    }
}
