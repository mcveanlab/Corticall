package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import java.util.Random;

public class LargeDelGenerator implements VariantGenerator {
    private int seqIndex;

    public LargeDelGenerator(int seqIndex) { this.seqIndex = seqIndex; }

    @Override
    public String getType() { return "LARGE_DEL"; }

    @Override
    public int getSeqIndex() { return seqIndex; }

    @Override
    public GeneratedVariant permute(String seq, int posIndex, Random rng) {
        int end = posIndex + rng.nextInt(1000) + 1000;

        String oldAllele = seq.substring(posIndex, posIndex + rng.nextInt(1000) + 1000);

        return new GeneratedVariant(getType(), getSeqIndex(), posIndex, oldAllele, "");
    }
}
