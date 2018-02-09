package uk.ac.ox.well.cortexjdk.utils.assembler;

import kotlin.sequences.Sequence;
import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinksIterable;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinksRecord;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

public class TempLinksAssemblerTest {
    private File findMcCortexBinary(int k) {
        File ctx = null;

        if (k > 0 && k <= 31) {
            ctx = new File("/Users/kiran/repositories/mccortex/bin/mccortex31");
        } else if (k > 31 && k <= 63) {
            ctx = new File("/Users/kiran/repositories/mccortex/bin/mccortex63");
        } else if (k > 63 && k <= 95) {
            ctx = new File("/Users/kiran/repositories/mccortex/bin/mccortex95");
        }

        return ctx != null && ctx.exists() ? ctx : null;
    }

    @Test
    public void testArbitraryLinksConstruction() {
        File ctxBin = findMcCortexBinary(31);

        if (ctxBin != null) {
            Random rng = new Random(0);

            List<String> sequences = new ArrayList<>();
//            for (int i = 0; i < 100; i++) {
//                String seq = new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(rng.nextInt(190) + 10));
//                sequences.add(seq);
//            }

            String repeat = new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(5));

            sequences.add(
                    new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(20)) +
                    repeat + repeat + repeat + repeat +
                    new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(20))
            );

            for (String seq : sequences) {
                Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
                haplotypes.put("test", Collections.singletonList(seq));

                CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 5);
                CortexLinks l = TempLinksAssembler.buildLinks(g, haplotypes, "test");

                try {
                    File tempFasta = File.createTempFile("tempseq", ".fa");
                    tempFasta.deleteOnExit();

                    File tempCtx = File.createTempFile("tempseq", ".ctx");
                    tempCtx.deleteOnExit();

                    File tempCtp = File.createTempFile("tempseq", ".ctp.gz");
                    tempCtp.deleteOnExit();

                    PrintWriter pw = new PrintWriter(tempFasta);
                    pw.println(">1");
                    pw.println(haplotypes.get("test").iterator().next());
                    pw.close();

                    Process p1 = Runtime.getRuntime().exec(ctxBin.getAbsolutePath() + " build -f -m 1G -k 5 -S -s test -1 " + tempFasta.getAbsolutePath() + " " + tempCtx.getAbsolutePath());
                    p1.waitFor();

                    Process p2 = Runtime.getRuntime().exec(ctxBin.getAbsolutePath() + " thread -f -m 1G -1 " + tempFasta.getAbsolutePath() + " -o " + tempCtp.getAbsolutePath() + " " + tempCtx.getAbsolutePath());
                    p2.waitFor();

                    CortexLinksIterable ll = new CortexLinksIterable(tempCtp);

                    int numLinks = 0;
                    for (CortexLinksRecord clr : ll) {
                        Assert.assertTrue(l.containsKey(clr.getKmer()));

                        numLinks += clr.getJunctions().size();
                    }

                    Assert.assertEquals(ll.getNumKmersWithLinks(), l.size());
                    Assert.assertEquals(ll.getNumLinks(), numLinks);
                } catch (IOException e) {
                    throw new CortexJDKException("Unable to create temp files", e);
                } catch (InterruptedException e) {
                    throw new CortexJDKException("Interrupted", e);
                }
            }
        }
    }

}
