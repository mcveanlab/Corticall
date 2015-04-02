package uk.ac.ox.well.indiana.commands.simulate;

import com.google.common.base.Joiner;
import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class FindVariantsInContigsTest {
    private final String KNOWN_TABLE = "testdata/child_0.reads_perfect.foundVariants.txt";
    private final String CONTIGS = "testdata/child_0.reads_perfect.contigs.metrics";

    private Map<String, Set<Map<String, String>>> loadExistingEvents() {
        Map<String, Set<Map<String, String>>> knownVariants = new HashMap<String, Set<Map<String, String>>>();

        TableReader tr = new TableReader(KNOWN_TABLE);

        for (Map<String, String> te : tr) {
            String contig1 = te.get("contig1");

            if (!knownVariants.containsKey(contig1)) {
                knownVariants.put(contig1, new HashSet<Map<String, String>>());
            }

            knownVariants.get(contig1).add(te);
        }

        return knownVariants;
    }

    private Map<String, String> loadContigs() {
        Map<String, String> contigs = new HashMap<String, String>();
        TableReader tr = new TableReader(CONTIGS);

        for (Map<String, String> te : tr) {
            String contigName = te.get("contigName");
            String seq = te.get("seq");

            contigs.put(contigName, seq);
        }

        return contigs;
    }

    @Test
    public void testAllelePositions() {
        Map<String, Set<Map<String, String>>> knownEvents = loadExistingEvents();
        Map<String, String> contigs = loadContigs();

        for (String contigName : knownEvents.keySet()) {
            for (Map<String, String> ke : knownEvents.get(contigName)) {
                if (ke.get("variantFound").equals("TRUE")) {
                    String contig1Name = ke.get("contig1");
                    String flank1 = ke.get("flank1");
                    String flank2 = ke.get("flank2");

                    int pos1 = Integer.valueOf(ke.get("pos1"));
                    int pos2 = Integer.valueOf(ke.get("pos2"));

                    String seq = contigs.get(contig1Name);

                    String testFlank1 = seq.substring(pos1, pos1 + flank1.length());
                    String testFlank2 = seq.substring(pos2, pos2 + flank2.length());

                    if (flank1.equals(testFlank1) && flank2.equals(testFlank2)) {
                        int variantPos = Integer.valueOf(ke.get("variantPos"));

                        Assert.assertEquals(variantPos, pos1 + flank1.length() + 1);
                        Assert.assertEquals(variantPos + ke.get("altAllele").length() - 1, pos2);
                    }
                }
            }
        }
    }
}
