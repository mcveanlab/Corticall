package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import java.util.Random;

public class SmallDelGenerator implements VariantGenerator {
    private int seqIndex;

    public SmallDelGenerator(int seqIndex) { this.seqIndex = seqIndex; }

    @Override
    public String getType() { return "SMALL_DEL"; }

    @Override
    public int getSeqIndex() { return seqIndex; }

    @Override
    public GeneratedVariant permute(String seq, int posIndex, Random rng) {
        String oldAllele = seq.substring(posIndex, posIndex + rng.nextInt(19) + 1);

        return new GeneratedVariant(getType(), getSeqIndex(), posIndex, oldAllele, "");
    }
}
