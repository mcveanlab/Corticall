package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.util.Random;

public class InvGenerator implements VariantGenerator {
    private int seqIndex;

    public InvGenerator(int seqIndex) { this.seqIndex = seqIndex; }

    @Override
    public String getType() { return "INV"; }

    @Override
    public int getSeqIndex() { return seqIndex; }

    @Override
    public GeneratedVariant permute(String seq, int posIndex, Random rng, int length) {
        String oldAllele = seq.substring(posIndex, posIndex + length);
        String newAllele = SequenceUtils.reverseComplement(oldAllele);

        return new GeneratedVariant(getType(), getSeqIndex(), posIndex, oldAllele, newAllele);
    }
}
