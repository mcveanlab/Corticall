package uk.ac.ox.well.cortexjdk.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.util.MathUtils;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.statistics.distributions.EmpiricalDistribution;

import java.util.*;

/**
 * Created by kiran on 19/11/2017.
 */
public class SimulateChild extends Module {
    @Argument(fullName="ref1", shortName="r1", doc="Ref 1")
    public IndexedFastaSequenceFile REF1;

    @Argument(fullName="ref2", shortName="r2", doc="Ref 2")
    public IndexedFastaSequenceFile REF2;

    @Argument(fullName="chrs1", shortName="c1", doc="Chrs 1")
    public ArrayList<String> CHRS1;

    @Argument(fullName="chrs2", shortName="c2", doc="Chrs 2")
    public ArrayList<String> CHRS2;

    @Argument(fullName="mu", shortName="m", doc="Mu array")
    public ArrayList<Double> MUS;

    private Random rng = new Random(System.currentTimeMillis());

    @Override
    public void execute() {
        if (CHRS1.size() != CHRS2.size() || CHRS1.size() != MUS.size()) {
            throw new CortexJDKException("Array size mismatch");
        }

        for (int i = 0; i < MUS.size(); i++) {
            String chrName1 = CHRS1.get(i);
            String chrName2 = CHRS2.get(i);

            double mu = MUS.get(i);
            double[] dp = poisson(mu, 10);

            EmpiricalDistribution ed = new EmpiricalDistribution(dp);

            String seq1 = REF1.getSequence(chrName1).getBaseString();
            String seq2 = REF2.getSequence(chrName2).getBaseString();

            Pair<List<String>, List<Integer>> res = recombine(seq1, seq2, ed.draw());
            List<String> pieces = res.getFirst();
            List<Integer> switches = res.getSecond();

            log.info("{} {} {} {} {} {}", i, mu, seq1.length(), seq2.length(), pieces.size(), Joiner.on("").join(switches));
        }
    }

    private Pair<List<String>, List<Integer>> addSNVs(List<String> seqs, int num) {
        List<StringBuilder> sbs = new ArrayList<>();
        for (String seq : seqs) {
            sbs.add(new StringBuilder(seq));
        }

        for (int i = 0; i < num; i++) {

        }

        return null;
    }

    private double[] poisson(double mu, int cmax) {
        double[] d = new double[cmax];

        for (int c = 0; c < cmax; c++) {
            d[c] = new PoissonDistribution(mu).probability(c);
        }

        return d;
    }

    private Pair<List<String>, List<Integer>> recombine(String seq1, String seq2, int numRecombs) {
        boolean order = rng.nextBoolean();
        String[] seqs = order ? new String[] { seq1, seq2 } : new String[] { seq2, seq1 };
        int[] parents = order ? new int[] { 1, 2 } : new int[] { 2, 1 };

        Map<String, Integer> kmersA = kmerize(seqs[0], 47);
        Map<String, Integer> kmersB = kmerize(seqs[1], 47);

        Map<String, Pair<Integer, Integer>> kmersO = overlap(kmersA, kmersB, 0.10f);
        List<String> kmers = new ArrayList<>(kmersO.keySet());

        List<String> pieces = new ArrayList<>();
        List<Integer> switches = new ArrayList<>();

        List<Pair<Integer, Integer>> recombPositions = new ArrayList<>();
        recombPositions.add(new Pair<>(0, 0));

        for (int r = 0; r < numRecombs; r++) {
            int index = rng.nextInt(kmers.size());

            recombPositions.add(kmersO.get(kmers.get(index)));
        }

        recombPositions.add(new Pair<>(seqs[0].length() - 1, seqs[1].length() - 1));
        recombPositions.sort(Comparator.comparing(Pair::getFirst));

        boolean drawFromFirst = true;
        for (int i = 0; i < recombPositions.size() - 1; i++) {
            Pair<Integer, Integer> p0 = recombPositions.get(i);
            Pair<Integer, Integer> p1 = recombPositions.get(i+1);

            if (drawFromFirst) {
                pieces.add(seqs[0].substring(p0.getFirst(), p1.getFirst()));
                switches.add(parents[0]);
            } else {
                pieces.add(seqs[1].substring(p0.getSecond(), p1.getSecond()));
                switches.add(parents[1]);
            }

            drawFromFirst = !drawFromFirst;
        }

        return new Pair<>(pieces, switches);
    }

    private Map<String, Integer> kmerize(String seq, int kmerSize) {
        Map<String, Integer> kmers = new LinkedHashMap<>();

        for (int i = 0; i <= seq.length() - kmerSize; i++) {
            String sk = seq.substring(i, i + kmerSize).toUpperCase();

            kmers.put(sk, i);
        }

        return kmers;
    }

    private Map<String, Pair<Integer, Integer>> overlap(Map<String, Integer> kmersA, Map<String, Integer> kmersB, float windowPct) {
        Set<String> kmers = new HashSet<>();
        kmers.addAll(kmersA.keySet());
        kmers.addAll(kmersB.keySet());

        Map<String, Pair<Integer, Integer>> kmersO = new LinkedHashMap<>();
        for (String kmer : kmers) {
            if (kmersA.containsKey(kmer) && kmersB.containsKey(kmer) && Math.abs(kmersA.get(kmer) - kmersB.get(kmer)) < (int) (windowPct*((float) kmersA.get(kmer)))) {
                kmersO.put(kmer, new Pair<>(kmersA.get(kmer), kmersB.get(kmer)));
            }
        }

        return kmersO;
    }
}
