package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import java.util.Random;

public interface VariantGenerator {
    String getType();
    int getSeqIndex();
    GeneratedVariant permute(String seq, int posIndex, Random rng, int length);
}

