package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.util.Random;

public class SmallInvGenerator implements VariantGenerator {
    private int seqIndex;

    public SmallInvGenerator(int seqIndex) { this.seqIndex = seqIndex; }

    @Override
    public String getType() { return "SMALL_INV"; }

    @Override
    public int getSeqIndex() { return seqIndex; }

    @Override
    public GeneratedVariant permute(String seq, int posIndex, Random rng) {
        int length = rng.nextInt(18) + 2;

        String oldAllele = seq.substring(posIndex, posIndex + length);
        String newAllele = SequenceUtils.reverseComplement(oldAllele);

        return new GeneratedVariant(getType(), getSeqIndex(), posIndex, oldAllele, newAllele);
    }
}
