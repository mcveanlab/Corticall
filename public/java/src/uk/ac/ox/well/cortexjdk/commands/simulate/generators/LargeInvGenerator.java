package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.util.Random;

public class LargeInvGenerator implements VariantGenerator {
    private int seqIndex;

    public LargeInvGenerator(int seqIndex) { this.seqIndex = seqIndex; }

    @Override
    public String getType() { return "LARGE_INV"; }

    @Override
    public int getSeqIndex() { return seqIndex; }

    @Override
    public GeneratedVariant permute(String seq, int posIndex, Random rng) {
        int end = posIndex + rng.nextInt(1000) + 1000;

        String oldAllele = seq.substring(posIndex, end);
        String newAllele = SequenceUtils.reverseComplement(oldAllele);

        return new GeneratedVariant(getType(), getSeqIndex(), posIndex, oldAllele, newAllele);
    }
}
