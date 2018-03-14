package uk.ac.ox.well.cortexjdk.utils.alignment.mosaic;

import htsjdk.samtools.util.StringUtil;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.util.Pair;
import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.util.*;

public class TestMosaicAligner {
    @Test
    public void testMosaicAlignment() {
        String[] templates = { new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(1000)),
                               new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(1000)) };

        List<Pair<String, String>> ks = new ArrayList<>();

        StringBuilder rb = new StringBuilder();
        int lastIndex = 0;
        int phase = 0;
        for (int recomb : Arrays.asList(200, 400, 600, 800)) {
            String subsection = templates[phase].substring(lastIndex, recomb);
            rb.append(subsection);

            ks.add(Pair.create("template" + phase, StringUtil.repeatCharNTimes(' ', lastIndex) + subsection));

            phase = phase == 0 ? 1 : 0;
            lastIndex = recomb;
        }
        ks.add(Pair.create("template" + phase, StringUtil.repeatCharNTimes(' ', 800) + templates[phase].substring(lastIndex, templates[0].length() - 1)));

        String subsection = templates[phase].substring(lastIndex, templates[0].length() - 1);
        rb.append(subsection);
        ks.add(0, Pair.create("query", rb.toString()));

        String query = rb.toString();
        Map<String, String> targets = new HashMap<>();
        targets.put("template0", templates[0]);
        targets.put("template1", templates[1]);

        MosaicAligner ma = new MosaicAligner();
        List<Triple<String, Pair<Integer, Integer>, String>> ps = ma.align(query, targets);

        Assert.assertEquals(ks.size(), ps.size());

        for (int i = 0; i < ks.size(); i++) {
            Assert.assertEquals(ps.get(i).getLeft(), ks.get(i).getFirst());
            Assert.assertTrue(unsharedKmers(ps.get(i).getRight().replaceAll(" ", ""), ks.get(i).getSecond().replaceAll(" ", ""), 47) <= 2);
        }
    }

    private int unsharedKmers(String s1, String s2, int k) {
        Set<String> ks1 = new HashSet<>();
        Set<String> ks2 = new HashSet<>();
        Set<String> u = new HashSet<>();

        for (int i = 0; i <= s1.length() - k; i++) {
            ks1.add(s1.substring(i, i + k));
            u.add(s1.substring(i, i + k));
        }

        for (int i = 0; i <= s2.length() - k; i++) {
            ks2.add(s2.substring(i, i + k));
            u.add(s2.substring(i, i + k));
        }

        int overlap = 0;
        for (String sk : u) {
            overlap += (ks1.contains(sk) && ks2.contains(sk)) ? 1 : 0;
        }

        return u.size() - overlap;
    }
}
