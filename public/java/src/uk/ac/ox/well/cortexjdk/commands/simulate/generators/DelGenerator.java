package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import java.util.Random;

public class DelGenerator implements VariantGenerator {
    private int seqIndex;

    public DelGenerator(int seqIndex) { this.seqIndex = seqIndex; }

    @Override
    public String getType() { return "DEL"; }

    @Override
    public int getSeqIndex() { return seqIndex; }

    @Override
    public GeneratedVariant permute(String seq, int posIndex, Random rng, int length) {
        String oldAllele = seq.substring(posIndex, posIndex + length + 1);
        String newAllele = seq.substring(posIndex, posIndex + 1);

        return new GeneratedVariant(getType(), getSeqIndex(), posIndex, oldAllele, newAllele);
    }
}
