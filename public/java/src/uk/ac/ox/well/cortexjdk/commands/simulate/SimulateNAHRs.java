package uk.ac.ox.well.cortexjdk.commands.simulate;

import com.google.common.base.Joiner;
import ngs.ReferenceSequence;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.gff.GFF3;
import uk.ac.ox.well.cortexjdk.utils.io.gff.GFF3Record;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by kiran on 26/10/2017.
 */
public class SimulateNAHRs extends Module {
    @Argument(fullName="ref1", shortName="r1", doc="Reference 1")
    public KmerLookup REF1;

    @Argument(fullName="gff1", shortName="g1", doc="GFF 1")
    public GFF3 GFF1;

    @Argument(fullName="ids1", shortName="i1", doc="IDs 1")
    public HashSet<String> IDS1;

    @Argument(fullName="ref2", shortName="r2", doc="Reference 2")
    public KmerLookup REF2;

    @Argument(fullName="gff2", shortName="g2", doc="GFF 2")
    public GFF3 GFF2;

    @Argument(fullName="ids2", shortName="i2", doc="IDs 2")
    public HashSet<String> IDS2;

    @Argument(fullName="seed", shortName="s", doc="Seed")
    public Long SEED = System.currentTimeMillis();

    @Argument(fullName="num", shortName="n", doc="Number of events to simulate")
    public Integer NUM = 1;

    @Output
    public PrintStream out;

    @Output(fullName="mout", shortName="mo", doc="M out")
    public PrintStream mout;

    private Random rng = new Random(SEED);

    @Override
    public void execute() {
        List<Pair<String, String>> vars = new ArrayList<>();
        vars.addAll(loadVars(REF1, GFF1, IDS1));
        vars.addAll(loadVars(REF2, GFF2, IDS2));

        for (int i = 0; i < NUM; i++) {
            int idx1 = 0, idx2 = 0;
            while (idx1 == idx2) {
                idx1 = rng.nextInt(vars.size());
                idx2 = rng.nextInt(vars.size());
            }

            Pair<String, Set<Integer>> newVar = recombine(vars, idx1, idx2);
            Pair<String, Set<Integer>> newVarMod = addSNVs(newVar.getFirst());

            out.println(">var_orig" + i + "\tvar1=" + vars.get(idx1).getFirst() + "_var2=" + vars.get(idx2).getFirst() + "_bks=" + Joiner.on(",").join(newVar.getSecond()));
            out.println(newVar.getFirst());

            mout.println(">var_mod" + i + "\tvar1=" + vars.get(idx1).getFirst() + "_var2=" + vars.get(idx2).getFirst() + "_bks=" + Joiner.on(",").join(newVar.getSecond()) + "_snvs=" + Joiner.on(",").join(newVarMod.getSecond()));
            mout.println(newVarMod.getFirst());
        }
    }

    private Pair<String, Set<Integer>> addSNVs(String seq) {
        StringBuilder sb = new StringBuilder(seq);
        int numSNVs = rng.nextInt(5);

        Set<Integer> pos = new TreeSet<>();
        for (int i = 0; i < numSNVs; i++) {
            int p = rng.nextInt(seq.length());

            char b;
            do {
                int bi = rng.nextInt(4);
                switch (bi) {
                    case 0: b = 'A'; break;
                    case 1: b = 'C'; break;
                    case 2: b = 'G'; break;
                    case 3: b = 'T'; break;
                    default: b = 'A'; break;
                }
            } while (seq.charAt(p) == b);

            sb.setCharAt(p, b);
            pos.add(p);
        }

        return new Pair<>(sb.toString(), pos);
    }

    private Pair<String, Set<Integer>> recombine(List<Pair<String, String>> vars, int idx1, int idx2) {
        Pair<String, String> vp1 = vars.get(idx1);
        Pair<String, String> vp2 = vars.get(idx2);

        String var1 = vp1.getSecond();
        String var2 = vp2.getSecond();

        String[] seqs = {var1, var2};
        int p = 0;
        int s = 0;

        int numRecombs = rng.nextInt(5) + 1;
        Set<Integer> recombPos = new TreeSet<>();
        for (int i = 0; i < numRecombs; i++) {
            int pos;
            do {
                pos = rng.nextInt(var1.length() < var2.length() ? var1.length() : var2.length());
            } while (withinWindow(recombPos, pos, 300));

            recombPos.add(pos);
            recombPos.add(pos + rng.nextInt(200) + 50);
        }

        StringBuilder sb = new StringBuilder();
        while (p < var1.length() && p < var2.length()) {
            if (s == 0) {
                sb.append(Character.toUpperCase(seqs[s].charAt(p)));
            } else {
                sb.append(Character.toLowerCase(seqs[s].charAt(p)));
            }

            p++;

            if (recombPos.contains(p)) {
                s = s == 0 ? 1 : 0;
            }
        }

        return new Pair<>(sb.toString(), recombPos);
    }

    private boolean withinWindow(Set<Integer> recombPos, int pos, int windowSize) {
        for (int i = pos - windowSize; i < pos + windowSize; i++) {
            if (recombPos.contains(i)) {
                return true;
            }
        }

        return false;
    }

    private List<Pair<String, String>> loadVars(KmerLookup kl, GFF3 gff, Set<String> ids) {
        Set<Pair<String, String>> vars = new HashSet<>();

        for (GFF3Record gr : gff) {
            if (ids.contains(gr.getAttribute("ID"))) {
                String seq = kl.getReferenceSequence().getSubsequenceAt(gr.getSeqid(), gr.getStart() - 500, gr.getEnd() + 500).getBaseString();

                if (gr.getStrand() == GFF3Record.Strand.NEGATIVE) {
                    seq = SequenceUtils.reverseComplement(seq);
                }

                vars.add(new Pair<>(gr.getAttribute("ID"), seq.toUpperCase()));
            }
        }

        return new ArrayList<>(vars);
    }
}
