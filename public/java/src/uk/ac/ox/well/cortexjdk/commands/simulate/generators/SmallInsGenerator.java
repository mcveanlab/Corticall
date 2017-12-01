package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.util.Random;

public class SmallInsGenerator implements VariantGenerator {
    private int seqIndex;

    public SmallInsGenerator(int seqIndex) { this.seqIndex = seqIndex; }

    @Override
    public String getType() { return "SMALL_INS"; }

    @Override
    public int getSeqIndex() { return seqIndex; }

    @Override
    public GeneratedVariant permute(String seq, int posIndex, Random rng) {
        String newAllele = new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(rng.nextInt(19) + 1));

        return new GeneratedVariant(getType(), getSeqIndex(), posIndex, seq.substring(posIndex, posIndex + 1), newAllele);
    }
}
