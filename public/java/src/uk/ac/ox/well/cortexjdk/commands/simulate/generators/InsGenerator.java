package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.util.Random;

public class InsGenerator implements VariantGenerator {
    private int seqIndex;

    public InsGenerator(int seqIndex) { this.seqIndex = seqIndex; }

    @Override
    public String getType() { return "INS"; }

    @Override
    public int getSeqIndex() { return seqIndex; }

    @Override
    public GeneratedVariant permute(String seq, int posIndex, Random rng, int length) {
        String oldAllele = seq.substring(posIndex, posIndex + 1);
        String newAllele = oldAllele + new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(length));

        return new GeneratedVariant(getType(), getSeqIndex(), posIndex, seq.substring(posIndex, posIndex + 1), newAllele);
    }
}
