package uk.ac.ox.well.cortexjdk.utils.alignment.graph;

import com.google.common.base.Joiner;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.StringUtil;
import org.apache.commons.math3.util.Pair;
import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.cortexjdk.Main;
import uk.ac.ox.well.cortexjdk.utils.assembler.TempGraphAssembler;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.util.*;

public class TestReadToGraphAligner {
    @Test
    public void testSingleEndAlignment() throws Exception {
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("mom", Arrays.asList("AGTTCGAATCTGGTCATGGTATTAGCGTATCGCCATCGATCAGGGCTATATGCT"));
        haplotypes.put("kid", Arrays.asList("AGTTCGAATCTGGTCATGGTCGATCAGGGCTATATGCT", "AGTTCGAATCTGGTCATGGTATTAGCGTATCGCCATCGATCAGGGCTATATGCT"));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 7);
        ReadToGraphAligner ga = new ReadToGraphAligner(g, null);

        Map<String, Boolean> reads = new LinkedHashMap<>();
        reads.put("AGTTCGAATCTGGTCATGGTCGATCAGGGCTATATGCT", true);
        reads.put("CGATCAGGGCTATATGCT", false);
        reads.put("AGTTCGAATCTG", false);
        reads.put("TCATGGTCGATCA", true);
        reads.put("ATGGTATTA", false);

        for (String readSequence : reads.keySet()) {
            FastqRecord fq = new FastqRecord("test", readSequence, null, StringUtil.repeatCharNTimes('I', readSequence.length()));

            SingleEndAlignmentInfo sai = ga.align(fq);
            boolean isSplitRead = sai.isSplitRead(1, 0);

            Assert.assertEquals(isSplitRead, (boolean) reads.get(readSequence));
        }
    }

    @Test
    public void testPairedEndAlignment() throws Exception {
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("mom", Arrays.asList(                                          "AGTTCGAATCTGGTCATGGTATTAGCGTATCGCCATCGATCAGGGCTATATGCT"));
        haplotypes.put("kid", Arrays.asList("AGTTCGAATCTGGTCATGGTCGATCAGGGCTATATGCT", "AGTTCGAATCTGGTCATGGTATTAGCGTATCGCCATCGATCAGGGCTATATGCT"));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 7);
        ReadToGraphAligner ga = new ReadToGraphAligner(g, null, 0, 1);

        Set<Pair<String, String>> reads = new LinkedHashSet<>();
        reads.add(new Pair<>("AGTTCGAATCTGG", SequenceUtils.reverseComplement("GGGCTATATGCT")));
        reads.add(new Pair<>("AGCATATAGCCCT", SequenceUtils.reverseComplement("CAGATTCGAACT")));

        for (Pair<String, String> read : reads) {
            FastqRecord fq1 = new FastqRecord("test", read.getFirst(), null, StringUtil.repeatCharNTimes('I', read.getFirst().length()));
            FastqRecord fq2 = new FastqRecord("test", read.getSecond(), null, StringUtil.repeatCharNTimes('I', read.getSecond().length()));

            PairedEndAlignmentInfo pai = ga.align(fq1, fq2);

            Main.getLogger().info("{} {} {}", read.getFirst(), read.getSecond(), Joiner.on(" ").withKeyValueSeparator(" => ").join(pai.getAlignmentDistanceMap()));
        }
    }

    @Test
    public void testPairedEndAlignmentWithSV() throws Exception {
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("mom", Arrays.asList("AGTTCGAATCTGGTCATGGTATTAGCGTATCGCCATCGATCAGGGCTATATGCT"));
        haplotypes.put("kid", Arrays.asList("AGTTCGAATCTGGTCATGGTATTAGGGCTAAATAGTTCATGATTGCATGTAGTCATGTGCGCATGACTAGTGCTCCCATAGCGTATCGCCATCGATCAGGGCTATATGCT"));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 11);
        ReadToGraphAligner ga = new ReadToGraphAligner(g, null, 0, 1);

        Set<Pair<String, String>> reads = new LinkedHashSet<>();
        reads.add(new Pair<>("AGTTCGAATCTGG", SequenceUtils.reverseComplement("GGGCTATATGCT")));

        for (Pair<String, String> read : reads) {
            FastqRecord fq1 = new FastqRecord("test", read.getFirst(), null, StringUtil.repeatCharNTimes('I', read.getFirst().length()));
            FastqRecord fq2 = new FastqRecord("test", read.getSecond(), null, StringUtil.repeatCharNTimes('I', read.getSecond().length()));

            PairedEndAlignmentInfo pai = ga.align(fq1, fq2);

            Main.getLogger().info("{} {} {}", read.getFirst(), SequenceUtils.reverseComplement(read.getSecond()), Joiner.on(" ").withKeyValueSeparator(" => ").join(pai.getAlignmentDistanceMap()));
        }
    }
}
