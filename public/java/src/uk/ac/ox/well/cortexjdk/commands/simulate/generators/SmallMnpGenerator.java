package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.util.Random;

public class SmallMnpGenerator implements VariantGenerator {
    private int seqIndex;

    public SmallMnpGenerator(int seqIndex) { this.seqIndex = seqIndex; }

    @Override
    public String getType() { return "SMALL_MNP"; }

    @Override
    public int getSeqIndex() { return seqIndex; }

    @Override
    public GeneratedVariant permute(String seq, int posIndex, Random rng) {
        int length = rng.nextInt(18) + 2;
        String newAllele = new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(length));

        return new GeneratedVariant(getType(), getSeqIndex(), posIndex, seq.substring(posIndex, posIndex + length), newAllele);
    }
}
