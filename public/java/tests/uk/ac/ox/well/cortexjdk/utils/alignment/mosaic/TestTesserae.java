package uk.ac.ox.well.cortexjdk.utils.alignment.mosaic;

import htsjdk.samtools.util.StringUtil;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.util.Pair;
import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.util.*;

public class TestTesserae {
    @Test
    public void testMosaicAlignment() {
        String[] templates = { new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(1000)),
                               new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(1000)) };

        List<Triple<String, String, Pair<Integer, Integer>>> ks = new ArrayList<>();

        StringBuilder rb = new StringBuilder();
        int lastIndex = 0;
        int phase = 0;
        for (int recomb : Arrays.asList(200, 400, 600, 800)) {
            String subsection = templates[phase].substring(lastIndex, recomb);
            rb.append(subsection);

            ks.add(Triple.of("template" + phase, StringUtil.repeatCharNTimes(' ', lastIndex) + subsection, Pair.create(lastIndex, recomb - 1)));

            phase = phase == 0 ? 1 : 0;
            lastIndex = recomb;
        }
        ks.add(Triple.of("template" + phase, StringUtil.repeatCharNTimes(' ', 800) + templates[phase].substring(lastIndex, templates[0].length() - 1), Pair.create(lastIndex, 999)));

        String subsection = templates[phase].substring(lastIndex, templates[0].length() - 1);
        rb.append(subsection);
        ks.add(0, Triple.of("query", rb.toString(), Pair.create(0, rb.length())));

        String query = rb.toString();
        Map<String, String> targets = new HashMap<>();
        targets.put("template0", templates[0]);
        targets.put("template1", templates[1]);

        Tesserae ma = new Tesserae();
        List<Triple<String, String, Pair<Integer, Integer>>> ps = ma.align(query, targets);

        Assert.assertEquals(ks.size(), ps.size());

        for (int i = 0; i < ks.size(); i++) {
            Assert.assertEquals(ps.get(i).getLeft(), ks.get(i).getLeft());
            Assert.assertTrue(unsharedKmers(ps.get(i).getMiddle().replaceAll(" ", ""), ks.get(i).getMiddle().replaceAll(" ", ""), 47) <= 2);

            System.out.println(ma);
        }
    }

    @Test
    private void smallTest() {
        String[] templates = { "GTAGGCGAGTCCCGTTTATA", "CCACAGAAGATGACGCCATT" };

        Map<String, String> targets = new HashMap<>();
        targets.put("template0", templates[0]);
        targets.put("template1", templates[1]);

        String query = "GTAGGCGAGATGACGCCAT";

        Tesserae ma = new Tesserae();
        List<Triple<String, String, Pair<Integer, Integer>>> ps = ma.align(query, targets);

        Assert.assertEquals(3, ps.size());

        /*
            query (0-18) GTAGGCGAGATGACGCCAT
                         |||||||||||||||||||
         template0 (0-6) GTAGGCG
        template1 (7-18)        AGATGACGCCAT
         */
    }

    @Test
    private void anotherSmallTest() {
        String query = "CGAACAGGATGTAGGCGAGATGACGCCATTTATTCTTTTCGTGCATAACAAAACGATAGTAG";

        String[] templates = {
                       "CGAACAGGATCAGGGATAAAACAAATTGATTATTCTTTTCGTGCATAACACGATAGTAG",
                       "GTCATACGACCGTAGGCGAGATGACGCCATTTATTACGGATATTATATTTATATA"
        };

        Map<String, String> targets = new HashMap<>();
        targets.put("template0", templates[0]);
        targets.put("template1", templates[1]);
        //targets.put("template2", templates[2]);
        //targets.put("template3", templates[3]);

        Tesserae ma = new Tesserae();
        List<Triple<String, String, Pair<Integer, Integer>>> ps = ma.align(query, targets);

        //Assert.assertEquals(3, ps.size());

        System.out.println(ma);

        /*
            query (0-18) GTAGGCGAGATGACGCCAT
                         |||||||||||||||||||
         template0 (0-6) GTAGGCG
        template1 (7-18)        AGATGACGCCAT
         */
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
