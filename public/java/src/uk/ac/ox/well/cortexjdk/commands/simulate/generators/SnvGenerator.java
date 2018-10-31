package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import java.util.Random;

public class SnvGenerator implements VariantGenerator {
    private int seqIndex;

    public SnvGenerator(int seqIndex) { this.seqIndex = seqIndex; }

    @Override
    public String getType() { return "SNV"; }

    @Override
    public int getSeqIndex() {
        return seqIndex;
    }

    @Override
    public GeneratedVariant permute(String seq, int posIndex, Random rng, int length) {
        final char[] bases = { 'A', 'C', 'G', 'T' };

        char base;
        do {
            base = bases[rng.nextInt(bases.length)];
        } while (base == Character.toUpperCase(seq.charAt(posIndex)) || base == Character.toLowerCase(seq.charAt(posIndex)));

        return new GeneratedVariant(getType(), getSeqIndex(), posIndex, seq.substring(posIndex, posIndex + 1), String.valueOf(base));
    }
}
