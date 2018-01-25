package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import java.util.Random;

public class TandemDuplicationGenerator implements VariantGenerator {
    private int seqIndex;

    public TandemDuplicationGenerator(int seqIndex) { this.seqIndex = seqIndex; }

    @Override
    public String getType() { return "TD"; }

    @Override
    public int getSeqIndex() { return seqIndex; }

    @Override
    public GeneratedVariant permute(String seq, int posIndex, Random rng) {
        int length = rng.nextInt(18) + 2;

        String oldAllele = seq.substring(posIndex, posIndex + length);
        String newAllele = oldAllele + oldAllele;

        return new GeneratedVariant(getType(), getSeqIndex(), posIndex, oldAllele, newAllele);
    }
}
