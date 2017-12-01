package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import org.apache.ivy.util.StringUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class StrExpGenerator implements VariantGenerator {
    private int seqIndex;

    public StrExpGenerator(int seqIndex) { this.seqIndex = seqIndex; }

    @Override
    public String getType() { return "STR_EXP"; }

    @Override
    public int getSeqIndex() { return seqIndex; }

    @Override
    public GeneratedVariant permute(String seq, int posIndex, Random rng) {
        List<Integer> loci = new ArrayList<>();

        int s = rng.nextInt(4) + 2;

        for (int i = 0; i < seq.length() - s; i++) {
            String repeatUnit = seq.substring(i, i + s);
            if (i + s + s < seq.length() && repeatUnit.equals(seq.substring(i + s, i + s + s))) {
                loci.add(i);
            }
        }

        int p = rng.nextInt(loci.size());
        int l = loci.get(p);

        String repeatUnit = seq.substring(l, l + s);
        int n = rng.nextInt(4) + 2;

        return new GeneratedVariant(getType(), getSeqIndex(), l, repeatUnit, StringUtils.repeat(repeatUnit, n));
    }
}
