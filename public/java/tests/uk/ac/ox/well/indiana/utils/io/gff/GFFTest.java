package uk.ac.ox.well.indiana.utils.io.gff;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class GFFTest {
    @Test
    public void loadGFFTest() {
        GFF3 gff = new GFF3("testdata/Pf3D7_14_v3.gff");

        List<GFF3Record> l1 = new ArrayList<>();
        List<GFF3Record> l2 = new ArrayList<>();

        for (GFF3Record gr : gff) {
            l1.add(gr);
        }

        for (GFF3Record gr : gff) {
            l2.add(gr);
        }

        Assert.assertEquals(l2.size(), l1.size());
    }
}
