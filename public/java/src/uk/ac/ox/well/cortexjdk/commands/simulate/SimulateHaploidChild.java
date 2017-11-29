package uk.ac.ox.well.cortexjdk.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.gff.GFF3;
import uk.ac.ox.well.cortexjdk.utils.io.gff.GFF3Record;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
import uk.ac.ox.well.cortexjdk.utils.statistics.distributions.EmpiricalDistribution;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by kiran on 19/11/2017.
 */
public class SimulateHaploidChild extends Module {
    @Argument(fullName="ref1", shortName="r1", doc="Ref 1")
    public IndexedFastaSequenceFile REF1;

    @Argument(fullName="ref2", shortName="r2", doc="Ref 2")
    public IndexedFastaSequenceFile REF2;

    @Argument(fullName="chrs1", shortName="c1", doc="Chrs 1")
    public ArrayList<String> CHRS1;

    @Argument(fullName="chrs2", shortName="c2", doc="Chrs 2")
    public ArrayList<String> CHRS2;

    @Argument(fullName="gff1", shortName="g1", doc="GFF 1")
    public GFF3 GFF1;

    @Argument(fullName="gff2", shortName="g2", doc="GFF 2")
    public GFF3 GFF2;

    @Argument(fullName="mu", shortName="m", doc="Mu array")
    public ArrayList<Double> MUS;

    @Argument(fullName="seed", shortName="s", doc="Random seed")
    public Long SEED = System.currentTimeMillis();

    @Argument(fullName="numBubbles", shortName="b", doc="Number of bubbles per type per chr")
    public Integer NUM_BUBBLES = 1;

    @Output
    public PrintStream out;

    @Output(fullName="variantsOut", shortName="vo", doc="Variants out")
    public PrintStream vout;

    private Random rng;

    private void initialize() {
        if (CHRS1.size() != CHRS2.size() || CHRS1.size() != MUS.size()) {
            throw new CortexJDKException("Array size mismatch");
        }
        rng = new Random(SEED);
    }

    @Override
    public void execute() {
        initialize();

        List<Pair<List<String>, List<Integer>>> seqs = new ArrayList<>();
        List<GFF3Record> gffs = new ArrayList<>();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing chromosomes")
                .message("processed")
                .maxRecord(CHRS1.size())
                .make(log);

        for (int i = 0; i < MUS.size(); i++) {
            String seq1 = REF1.getSequence(CHRS1.get(i)).getBaseString().toUpperCase();
            String seq2 = REF2.getSequence(CHRS2.get(i)).getBaseString().toLowerCase();

            double mu = MUS.get(i);
            double[] dp = poisson(mu, 10);
            EmpiricalDistribution ed = new EmpiricalDistribution(dp, rng);
            int numRecombs = ed.draw();

            Pair<List<String>, List<Integer>> res = recombine(seq1, seq2, numRecombs, 47, 0.10f);
            gffs.addAll(liftover(res.getFirst(), i, CHRS1.get(i), CHRS2.get(i)));

            seqs.add(res);

            pm.update();
        }

        Set<GeneratedVariant> vs = new TreeSet<>();

        vs.addAll(makeNAHR(seqs, gffs));

        for (int i = 0; i < CHRS1.size(); i++) {
            Pair<List<String>, List<Integer>> res = seqs.get(i);
            List<String> pieces = res.getFirst();

            vs.addAll(makeBubble(pieces, NUM_BUBBLES, rng, new SnvGenerator(i)));
            vs.addAll(makeBubble(pieces, NUM_BUBBLES, rng, new SmallInsGenerator(i)));
            vs.addAll(makeBubble(pieces, NUM_BUBBLES, rng, new SmallDelGenerator(i)));
            vs.addAll(makeBubble(pieces, NUM_BUBBLES, rng, new SmallMnpGenerator(i)));
        }

        List<String> newSeqs = collapse(seqs, vs);

        for (int i = 0; i < newSeqs.size(); i++) {
            out.println(">" + i);
            out.println(newSeqs.get(i));
        }
    }

    private List<String> collapse(List<Pair<List<String>, List<Integer>>> seqs, Set<GeneratedVariant> vs) {
        List<GeneratedVariant> lvs = new ArrayList<>(vs);

        List<StringBuilder> newSbs = new ArrayList<>();
        for (int i = 0; i < seqs.size(); i++) {
            StringBuilder sbs = new StringBuilder();
            for (int j = 0; j < seqs.get(i).getSecond().size(); j++) {
                String sb = seqs.get(i).getFirst().get(j);
                int sw = seqs.get(i).getSecond().get(j);

                sb = sw == 1 ? sb.toUpperCase() : sb.toLowerCase();

                sbs.append(sb);
            }

            newSbs.add(sbs);
        }

        vout.println(Joiner.on("\t").join("index", "chr", "pos", "type", "old", "new", "sleft", "sright"));
        for (int i = lvs.size() - 1; i >= 0; i--) {
            GeneratedVariant gv = lvs.get(i);
            StringBuilder sb = newSbs.get(gv.seqIndex);

            log.info("{} {} {} {}", i, sb.length(), sb.substring(gv.posIndex, gv.posIndex + gv.oldAllele.length()), gv);

            sb = sb.replace(gv.posIndex, gv.posIndex + gv.oldAllele.length(), gv.newAllele);

            String seedLeft = sb.substring(gv.posIndex - 100, gv.posIndex);
            String seedRight = sb.substring(gv.posIndex + gv.newAllele.length(), gv.posIndex + gv.newAllele.length() + 100);

            vout.println(Joiner.on("\t").join(i,
                                              gv.getSeqIndex(),
                                              gv.getPosIndex(),
                                              gv.getType(),
                                              gv.getOldAllele().length() == 0 ? "." : gv.getOldAllele(),
                                              gv.getNewAllele().length() == 0 ? "." : gv.getNewAllele(),
                                              seedLeft,
                                              seedRight));
        }

        List<String> newSeqs = new ArrayList<>();
        for (StringBuilder sb : newSbs) {
            newSeqs.add(sb.toString());
        }

        return newSeqs;
    }

    private Set<GeneratedVariant> makeNAHR(List<Pair<List<String>, List<Integer>>> seqs, List<GFF3Record> gffs) {
        Set<GeneratedVariant> vs = new TreeSet<>();

        int vi1 = 0, vi2 = 0;
        while (vi1 == vi2) {
            vi1 = rng.nextInt(gffs.size());
            vi2 = rng.nextInt(gffs.size());
        }

        GFF3Record g1 = gffs.get(vi1);
        String v1  = Joiner.on("").join(seqs.get(Integer.valueOf(g1.getSeqid())).getFirst()).substring(g1.getStart(), g1.getEnd());
        String vo1 = g1.getStrand() == GFF3Record.Strand.POSITIVE ? v1 : SequenceUtils.reverseComplement(v1);

        GFF3Record g2 = gffs.get(vi2);
        String v2  = Joiner.on("").join(seqs.get(Integer.valueOf(g2.getSeqid())).getFirst()).substring(g2.getStart(), g2.getEnd());
        String vo2 = g2.getStrand() == GFF3Record.Strand.POSITIVE ? v2 : SequenceUtils.reverseComplement(v2);

        int numVarRecombs = rng.nextInt(5) + 1;
        Pair<List<String>, List<Integer>> vres = recombine(vo1, vo2, numVarRecombs,21,.25f);

        vs.add(new GeneratedVariant("NAHR-INS", Integer.valueOf(g1.getSeqid()), g1.getStart(), v1, Joiner.on("").join(vres.getFirst())));
        vs.add(new GeneratedVariant("NAHR-DEL", Integer.valueOf(g2.getSeqid()), g2.getStart(), v2, ""));

        return vs;
    }

    private Set<GFF3Record> liftover(List<String> pieces, int seqIndex, String chr1, String chr2) {
        Set<GFF3Record> grs = new TreeSet<>();

        String seq = Joiner.on("").join(pieces);

        for (GFF3Record oldGr : GFF1.getContained(new Interval(chr1, 0, REF1.getSequence(chr1).length()))) {
            String gene = REF1.getSubsequenceAt(oldGr.getSeqid(), oldGr.getStart(), oldGr.getEnd()).getBaseString();

            int newIndex = seq.indexOf(gene);

            if (newIndex >= 0) {
                GFF3Record newGr = new GFF3Record(oldGr);
                newGr.setSeqid(String.valueOf(seqIndex));
                newGr.setStart(newIndex + 1);
                newGr.setEnd(newIndex + gene.length());

                grs.add(newGr);
            }
        }

        for (GFF3Record oldGr : GFF2.getContained(new Interval(chr2, 0, REF2.getSequence(chr2).length()))) {
            String gene = REF2.getSubsequenceAt(oldGr.getSeqid(), oldGr.getStart(), oldGr.getEnd()).getBaseString();

            int newIndex = seq.indexOf(gene);

            if (newIndex >= 0) {
                GFF3Record newGr = new GFF3Record(oldGr);
                newGr.setSeqid(String.valueOf(seqIndex));
                newGr.setStart(newIndex + 1);
                newGr.setEnd(newIndex + gene.length());

                grs.add(newGr);
            }
        }

        return grs;
    }

    private Set<GeneratedVariant> makeBubble(List<String> seqs, int num, Random rng, VariantGenerator v) {
        Set<GeneratedVariant> vs = new TreeSet<>();

        StringBuilder sb = new StringBuilder();
        for (String seq : seqs) {
            sb.append(seq);
        }
        String seq = sb.toString();

        for (int i = 0; i < num; i++) {
            int posIndex = rng.nextInt(seq.length());

            GeneratedVariant gv = v.permute(seq, posIndex, rng);
            vs.add(gv);
        }

        return vs;
    }

    private double[] poisson(double mu, int cmax) {
        double[] d = new double[cmax];

        for (int c = 0; c < cmax; c++) {
            d[c] = new PoissonDistribution(mu).probability(c);
        }

        return d;
    }

    private Pair<List<String>, List<Integer>> recombine(String seq1, String seq2, int numRecombs, int kmerSize, float windowPct) {
        boolean order = rng.nextBoolean();
        String[] seqs = order ? new String[] { seq1, seq2 } : new String[] { seq2, seq1 };
        int[] parents = order ? new int[] { 1, 2 } : new int[] { 2, 1 };

        Map<String, Integer> kmersA = kmerize(seqs[0], kmerSize);
        Map<String, Integer> kmersB = kmerize(seqs[1], kmerSize);

        Map<String, Pair<Integer, Integer>> kmersO = overlap(kmersA, kmersB, windowPct);
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
