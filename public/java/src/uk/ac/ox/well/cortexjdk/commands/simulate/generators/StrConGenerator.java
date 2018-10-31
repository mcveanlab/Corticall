package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class StrConGenerator implements VariantGenerator {
    private int seqIndex;

    public StrConGenerator(int seqIndex) { this.seqIndex = seqIndex; }

    @Override
    public String getType() { return "STR_CON"; }

    @Override
    public int getSeqIndex() { return seqIndex; }

    @Override
    public GeneratedVariant permute(String seq, int posIndex, Random rng, int length) {
        List<Integer> loci = new ArrayList<>();

        int s = rng.nextInt(4) + 2;

        for (int i = 0; i < seq.length() - s; i++) {
            String repeatUnit = seq.substring(i, i + s);
            if (!repeatUnit.contains("N") && i + s + s < seq.length() && repeatUnit.equals(seq.substring(i + s, i + s + s))) {
                loci.add(i);
            }
        }

        int p = rng.nextInt(loci.size());
        int l = loci.get(p);

        String repeatUnit = seq.substring(l, l + s);
        int numAdjacentRepeats = 0;

        for (int i = l; i < seq.length() - s; i += s) {
            String ru = seq.substring(i, i + s);
            if (repeatUnit.equals(ru)) {
                numAdjacentRepeats++;
            } else {
                break;
            }
        }

        for (int i = l - s; i >= 0; i -= s) {
            String ru = seq.substring(i, i + s);
            if (repeatUnit.equals(ru)) {
                numAdjacentRepeats++;
            } else {
                break;
            }
        }

        int n = rng.nextInt(numAdjacentRepeats - 1) + 2;

        String oldAllele = seq.substring(l, l + (n*s));
        String newAllele = seq.substring(l, l + s);

        return new GeneratedVariant(getType(), getSeqIndex(), l, oldAllele, newAllele);
    }
}
