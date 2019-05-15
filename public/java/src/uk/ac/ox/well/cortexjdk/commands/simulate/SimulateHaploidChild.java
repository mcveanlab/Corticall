package uk.ac.ox.well.cortexjdk.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.commands.simulate.generators.*;
import uk.ac.ox.well.cortexjdk.utils.alignment.sw.SmithWaterman;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.gff.GFF3;
import uk.ac.ox.well.cortexjdk.utils.io.gff.GFF3Record;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexBinaryKmer;
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
    @Argument(fullName="parents", shortName="p", doc="Parent names")
    public ArrayList<String> PARENTS;

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

    @Argument(fullName="numVariants", shortName="v", doc="Number of variants")
    public Integer NUM_VARIANTS = 3;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Output
    public PrintStream out;

    @Output(fullName="fastaOut1", shortName="fo1", doc="Fasta out 1")
    public PrintStream fout1;

    @Output(fullName="fastaOut2", shortName="fo2", doc="Fasta out 2")
    public PrintStream fout2;

    @Output(fullName="variantsOut", shortName="vo", doc="Variants out")
    public PrintStream vout;

    @Output(fullName="kmersOut", shortName="ko", doc="Kmers out")
    public PrintStream kout;

    private Random rng;

    private void initialize() {
        if (CHRS1.size() != CHRS2.size() || CHRS1.size() != MUS.size()) {
            throw new CortexJDKException("Array size mismatch");
        }
        rng = new Random(SEED);

        vout.println(Joiner.on("\t").join("index", "chr", "start", "stop", "parent", "type", "old", "new", "sleft", "sright", "refChr", "refStart", "refStop"));
    }

    @Override
    public void execute() {
        initialize();

        List<Pair<List<String>, List<Integer>>> seqs = new ArrayList<>();
        List<GFF3Record> gffs = new ArrayList<>();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing chromosomes...")
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

            Pair<List<String>, List<Integer>> res = recombine(seq1, seq2, numRecombs, KMER_SIZE, 0.10f);
            gffs.addAll(liftover(res.getFirst(), i, CHRS1.get(i), CHRS2.get(i)));

            seqs.add(res);

            int start = 0;
            for (int j = 0; j < res.getFirst().size(); j++) {
                String piece = res.getFirst().get(j);
                int sw = res.getSecond().get(j);

                vout.println(Joiner.on("\t").join(-1, i+1, start, start + piece.length(), PARENTS.get(sw - 1), "RECOMB", ".", ".", ".", ".", ".", ".", "."));

                start += piece.length();
            }

            pm.update();
        }

        log.info("Adding variants...");
        Set<GeneratedVariant> vs = new TreeSet<>();
        Set<Integer> seqIndices = new HashSet<>();

        //vs.addAll(makeNAHR(seqs, gffs));

        for (int numBubble = 0; numBubble < NUM_VARIANTS; numBubble++) {
            int chrIndex = rng.nextInt(CHRS1.size());
            Pair<List<String>, List<Integer>> res = seqs.get(chrIndex);
            List<String> pieces = res.getFirst();

            int variantType = rng.nextInt(10);
            Set<GeneratedVariant> v;

            switch (variantType) {
                case 0:
                    v = makeBubble(pieces, 1, rng, new InsGenerator(chrIndex));
                    break;
                case 1:
                    v = makeBubble(pieces, 1, rng, new StrExpGenerator(chrIndex));
                    break;
                case 2:
                    v = makeBubble(pieces, 1, rng, new TandemDuplicationGenerator(chrIndex));
                    break;
                case 3:
                    v = makeBubble(pieces, 1, rng, new DelGenerator(chrIndex));
                    break;
                case 4:
                    v = makeBubble(pieces, 1, rng, new StrConGenerator(chrIndex));
                    break;
                case 5:
                    v = makeBubble(pieces, 1, rng, new MnpGenerator(chrIndex));
                    break;
                case 6:
                    v = makeBubble(pieces, 1, rng, new InvGenerator(chrIndex));
                    break;
                case 7:
                    v = makeNAHR(seqs, gffs);
                    break;
                case 8:
                default:
                    v = makeBubble(pieces, 1, rng, new SnvGenerator(chrIndex));
            }

            vs.addAll(v);
        }

        log.info("Collapsing into new linear haploid reference...");
        Set<CortexBinaryKmer> cbks = getParentalKmers(REF1, REF2, KMER_SIZE);
        List<String> newSeqs = collapse(seqs, cbks, vs);

        for (int i = 0; i < newSeqs.size(); i++) {
            if (seqIndices.contains(i)) {
                out.println(">chr" + (i + 1));
                out.println(newSeqs.get(i));

                ReferenceSequence r1 = REF1.getSequence(REF1.getSequenceDictionary().getSequence(i).getSequenceName());
                ReferenceSequence r2 = REF2.getSequence(REF2.getSequenceDictionary().getSequence(i).getSequenceName());

                fout1.println(">" + r1.getName());
                fout1.println(r1.getBaseString());
                fout2.println(">" + r2.getName());
                fout2.println(r2.getBaseString());
            }
        }

    }

    private Set<CortexBinaryKmer> getParentalKmers(IndexedFastaSequenceFile ref1, IndexedFastaSequenceFile ref2, int kmerSize) {
        ProgressMeter pm = new ProgressMeterFactory()
                .header("Computing parental kmers...")
                .message("processed")
                .maxRecord(ref1.getSequenceDictionary().size() + ref2.getSequenceDictionary().size())
                .make(log);

        Set<CortexBinaryKmer> cbks = new HashSet<>();
        for (IndexedFastaSequenceFile ref : Arrays.asList(ref1, ref2)) {
            ReferenceSequence rseq;
            while ((rseq = ref.nextSequence()) != null) {
                byte[] bseq = rseq.getBases();

                for (int i = 0; i <= bseq.length - kmerSize; i++) {
                    byte[] bk = new byte[kmerSize];
                    System.arraycopy(bseq, i, bk, 0, kmerSize);

                    boolean hasNs = false;
                    for (int j = 0; j < bk.length; j++) {
                        if (bk[j] == 'N' || bk[j] == 'n') {
                            hasNs = true;
                            break;
                        }
                    }

                    if (!hasNs) {
                        cbks.add(new CortexBinaryKmer(bk));
                    }
                }

                pm.update();
            }
        }

        return cbks;
    }

    private List<String> collapse(List<Pair<List<String>, List<Integer>>> seqs, Set<CortexBinaryKmer> cbks, Set<GeneratedVariant> vs) {
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

        for (int i = lvs.size() - 1; i >= 0; i--) {
            GeneratedVariant gv = lvs.get(i);
            StringBuilder sb = newSbs.get(gv.seqIndex);

            //log.info("{} {} {}", i, gv, sb.substring(gv.posIndex, gv.posIndex + gv.oldAllele.length()));

            String oldAlleleCheck = sb.substring(gv.posIndex, gv.posIndex + gv.oldAllele.length());
            if (!oldAlleleCheck.equals(gv.oldAllele)) {
                throw new CortexJDKException("WTF?");
            }

            sb = sb.replace(gv.posIndex, gv.posIndex + gv.oldAllele.length(), gv.newAllele);

            String seedLeft = sb.substring(gv.posIndex - 100, gv.posIndex);
            String seedRight = sb.substring(gv.posIndex + gv.newAllele.length(), gv.posIndex + gv.newAllele.length() + 100);

            Set<String> novelKmers = new LinkedHashSet<>();
            for (int p = gv.posIndex - 100; p <= gv.posIndex + gv.newAllele.length() + 100 - KMER_SIZE; p++) {
                String sk = sb.substring(p, p + KMER_SIZE).toUpperCase();
                CortexBinaryKmer cbk = new CortexBinaryKmer(sk.getBytes());

                if (!cbks.contains(cbk)) {
                    novelKmers.add(sk);
                }
            }

            IndexedFastaSequenceFile parentRef = Character.isUpperCase(seedLeft.charAt(0)) ? REF1 : REF2;
            ReferenceSequence parentSeq = parentRef.getSequence(parentRef.getSequenceDictionary().getSequence(gv.getSeqIndex()+1).getSequenceName());
            int refPosLeft = parentSeq.getBaseString().toUpperCase().indexOf(seedLeft.toUpperCase()) + seedLeft.length();
            int refPosRight = parentSeq.getBaseString().toUpperCase().indexOf(seedRight.toUpperCase()) + 1;

            if (gv.getLoci() != null) {
                List<Pair<Interval, Interval>> loci = gv.getLoci();

                for (Pair<Interval, Interval> locus : loci) {
                    vout.println(Joiner.on("\t").join(
                            i,
                            gv.getSeqIndex() + 1,
                            gv.getPosIndex(),
                            gv.getPosIndex() + gv.newAllele.length(),
                            Character.isUpperCase(seedLeft.charAt(0)) ? PARENTS.get(0) : PARENTS.get(1),
                            gv.getType(),
                            gv.getOldAllele().length() == 0 ? "." : gv.getOldAllele(),
                            gv.getNewAllele().length() == 0 ? "." : gv.getNewAllele(),
                            seedLeft,
                            seedRight,
                            parentSeq.getName(),
                            locus.getFirst().getContig() + ":" + locus.getFirst().getStart() + "-" + locus.getFirst().getEnd() + ":" + (locus.getFirst().isNegativeStrand() ? "-" : "+"),
                            locus.getSecond().getContig() + ":" + locus.getSecond().getStart() + "-" + locus.getSecond().getEnd() + ":" + (locus.getSecond().isNegativeStrand() ? "-" : "+")
                    ));
                }
            } else {
                vout.println(Joiner.on("\t").join(
                        i,
                        gv.getSeqIndex() + 1,
                        gv.getPosIndex(),
                        gv.getPosIndex() + gv.newAllele.length(),
                        Character.isUpperCase(seedLeft.charAt(0)) ? PARENTS.get(0) : PARENTS.get(1),
                        gv.getType(),
                        gv.getOldAllele().length() == 0 ? "." : gv.getOldAllele(),
                        gv.getNewAllele().length() == 0 ? "." : gv.getNewAllele(),
                        seedLeft,
                        seedRight,
                        parentSeq.getName(),
                        refPosLeft,
                        refPosRight
                ));
            }

            int nki = 0;
            for (String nk : novelKmers) {
                kout.println(Joiner.on("\t").join(i, novelKmers.size(), nki, nk));
                nki++;
            }
        }

        List<String> newSeqs = new ArrayList<>();
        for (StringBuilder sb : newSbs) {
            newSeqs.add(sb.toString());
        }

        return newSeqs;
    }

    private Set<GeneratedVariant> makeNAHR(List<Pair<List<String>, List<Integer>>> seqs, List<GFF3Record> gffs) {
        Set<GeneratedVariant> vs = new TreeSet<>();
        Triple<List<String>, List<Integer>, List<Pair<Interval, Interval>>> vres = null;

        do {
            int vi1 = 0, vi2 = 0;
            while (vi1 == vi2) {
                vi1 = rng.nextInt(gffs.size());
                vi2 = rng.nextInt(gffs.size());
            }

            GFF3Record g1 = gffs.get(vi1);
            String v1 = Joiner.on("").join(seqs.get(Integer.valueOf(g1.getSeqid())).getFirst()).substring(g1.getStart(), g1.getEnd());

            GFF3Record g2 = gffs.get(vi2);
            String v2 = Joiner.on("").join(seqs.get(Integer.valueOf(g2.getSeqid())).getFirst()).substring(g2.getStart(), g2.getEnd());

            if (g1.getInterval().length() > 500 && g2.getInterval().length() > 500) {
                int numVarRecombs = rng.nextInt(5) + 2;
                vres = recombine(seqs, g1, g2, numVarRecombs);

                List<Integer> switches = new ArrayList<>();
                switches.add(0);
                int sw = 0;
                for (String v : vres.getLeft()) {
                    sw += v.length();
                    switches.add(sw);
                }

                for (int i = 0; i < vres.getLeft().size(); i++) {
                    if (vres.getMiddle().get(i) == 1) {
                        log.info("{} seq1 {}", i, g1.getStart());
                    } else {
                        log.info("{} seq1 {}", i, g2.getStart());
                    }
                }

                String nalt = Joiner.on("").join(vres.getLeft());

                if (vres.getMiddle().get(0) == 1) {
                    if (g1.getStrand() == GFF3Record.Strand.NEGATIVE) {
                        nalt = SequenceUtils.reverseComplement(nalt);
                    }

                    SmithWaterman smw = new SmithWaterman();
                    String[] a = smw.getAlignment(v1, nalt);
                    log.info("a0: {}", a[0]);
                    log.info("a1: {}", a[1]);

                    vs.add(new GeneratedVariant("NAHR-INS", Integer.valueOf(g1.getSeqid()), g1.getStart(), v1, nalt, vres.getRight()));
                    vs.add(new GeneratedVariant("NAHR-DEL", Integer.valueOf(g2.getSeqid()), g2.getStart(), v2, ""));
                } else {
                    if (g2.getStrand() == GFF3Record.Strand.NEGATIVE) {
                        nalt = SequenceUtils.reverseComplement(nalt);
                    }

                    SmithWaterman smw = new SmithWaterman();
                    String[] a = smw.getAlignment(v2, nalt);
                    log.info("a0: {}", a[0]);
                    log.info("a1: {}", a[1]);

                    vs.add(new GeneratedVariant("NAHR-INS", Integer.valueOf(g2.getSeqid()), g2.getStart(), v2, nalt, vres.getRight()));
                    vs.add(new GeneratedVariant("NAHR-DEL", Integer.valueOf(g1.getSeqid()), g1.getStart(), v1, ""));
                }
            }
        } while (vres == null);

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
            GeneratedVariant gv;
            do {
                int length = rng.nextInt(1000);
                int posIndex = rng.nextInt(seq.length() - 2*length) + length;
                gv = v.permute(seq, posIndex, rng, length);
            } while (gv.getOldAllele().contains("N") || gv.getNewAllele().contains("N"));

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

    private Triple<List<String>, List<Integer>, List<Pair<Interval, Interval>>> recombine(List<Pair<List<String>, List<Integer>>> seqs, GFF3Record g1, GFF3Record g2, int numRecombs) {
        String v1 = Joiner.on("").join(seqs.get(Integer.valueOf(g1.getSeqid())).getFirst()).substring(g1.getStart(), g1.getEnd());
        String v1ref = g1.getSeqid();
        int v1start = g1.getStart();
        int v1end = g1.getEnd();
        for (IndexedFastaSequenceFile ref : Arrays.asList(REF1, REF2)) {
            ReferenceSequence rseq = ref.getSequence(ref.getSequenceDictionary().getSequence(Integer.valueOf(g1.getSeqid())+1).getSequenceName());

            int index = rseq.getBaseString().indexOf(v1);
            if (index >= 0) {
                v1ref = rseq.getName();
                v1start = index;
                v1end = v1start + v1.length();
                break;
            }
        }

        Interval[] i1 = new Interval[v1.length() + 1];
        for (int i = 0; i < v1.length() + 1; i++) {
            if (g1.getStrand() == GFF3Record.Strand.POSITIVE) {
                i1[i] = new Interval(v1ref, v1start + i, v1start + i, false, ".");
            } else {
                i1[i] = new Interval(v1ref, v1end - i, v1end - i, true, ".");
            }
        }
        String seq1 = (g1.getStrand() == GFF3Record.Strand.POSITIVE) ? v1 : SequenceUtils.reverseComplement(v1);

        String v2 = Joiner.on("").join(seqs.get(Integer.valueOf(g2.getSeqid())).getFirst()).substring(g2.getStart(), g2.getEnd());
        String v2ref = g2.getSeqid();
        int v2start = g2.getStart();
        int v2end = g2.getEnd();
        for (IndexedFastaSequenceFile ref : Arrays.asList(REF1, REF2)) {
            ReferenceSequence rseq = ref.getSequence(ref.getSequenceDictionary().getSequence(Integer.valueOf(g2.getSeqid())+1).getSequenceName());

            int index = rseq.getBaseString().indexOf(v2);
            if (index >= 0) {
                v2ref = rseq.getName();
                v2start = index;
                v2end = v2start + v2.length();
                break;
            }
        }

        Interval[] i2 = new Interval[v1.length() + 1];
        for (int i = 0; i < v1.length() + 1; i++) {
            if (g2.getStrand() == GFF3Record.Strand.POSITIVE) {
                i2[i] = new Interval(v2ref, v2start + i, v2start + i, false, ".");
            } else {
                i2[i] = new Interval(v2ref, v2end - i, v2end - i, true, ".");
            }
        }
        String seq2 = (g2.getStrand() == GFF3Record.Strand.POSITIVE) ? v2 : SequenceUtils.reverseComplement(v2);

        boolean order = rng.nextBoolean();
        String[] ss = order ? new String[] { seq1, seq2 } : new String[] { seq2, seq1 };
        int[] parents = order ? new int[] { 1, 2 } : new int[] { 2, 1 };
        Interval[] ii1 = order ? i1 : i2;
        Interval[] ii2 = order ? i2 : i1;

        List<Pair<Integer, Integer>> recombPositions = new ArrayList<>();
        recombPositions.add(new Pair<>(0, 0));

        int minLength = Math.min(seq1.length(), seq2.length());
        int recombLength = (int) Math.floor(minLength / (numRecombs+1));

        for (int r = 0; r < numRecombs; r++) {
            Pair<Integer, Integer> p = new Pair<>((r+1)*recombLength, (r+1)*recombLength);

            recombPositions.add(p);
        }

        recombPositions.add(new Pair<>(ss[0].length(), ss[1].length()));
        recombPositions.sort(Comparator.comparing(Pair::getFirst));

        List<String> pieces = new ArrayList<>();
        List<Integer> switches = new ArrayList<>();
        List<Pair<Interval, Interval>> loci = new ArrayList<>();

        boolean drawFromFirst = true;
        for (int i = 0; i < recombPositions.size() - 1; i++) {
            Pair<Integer, Integer> p0 = recombPositions.get(i);
            Pair<Integer, Integer> p1 = recombPositions.get(i+1);
            Pair<Interval, Interval> l;

            if (drawFromFirst) {
                pieces.add(ss[0].substring(p0.getFirst(), p1.getFirst()));
                switches.add(parents[0]);
                l = new Pair<>(ii1[p1.getFirst()], ii2[p1.getFirst()]);
            } else {
                pieces.add(ss[1].substring(p0.getSecond(), p1.getSecond()));
                switches.add(parents[1]);
                l = new Pair<>(ii2[p1.getSecond()], ii1[p1.getSecond()]);
            }

            drawFromFirst = !drawFromFirst;

            if (i < recombPositions.size() - 2) {
                loci.add(l);
            }
        }

        return Triple.of(pieces, switches, loci);
    }

    private Pair<List<String>, List<Integer>> recombine(String seq1, String seq2, int numRecombs) {
        boolean order = rng.nextBoolean();
        String[] seqs = order ? new String[] { seq1, seq2 } : new String[] { seq2, seq1 };
        int[] parents = order ? new int[] { 1, 2 } : new int[] { 2, 1 };

        List<Pair<Integer, Integer>> recombPositions = new ArrayList<>();
        recombPositions.add(new Pair<>(0, 0));

        int minLength = Math.min(seq1.length(), seq2.length());
        int recombLength = (int) Math.floor(minLength / (numRecombs+1));

        for (int r = 0; r < numRecombs; r++) {
            Pair<Integer, Integer> p = new Pair<>((r+1)*recombLength, (r+1)*recombLength);

            recombPositions.add(p);
        }

        recombPositions.add(new Pair<>(seqs[0].length() - 1, seqs[1].length() - 1));
        recombPositions.sort(Comparator.comparing(Pair::getFirst));

        List<String> pieces = new ArrayList<>();
        List<Integer> switches = new ArrayList<>();

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

    private Pair<List<String>, List<Integer>> recombine(String seq1, String seq2, int numRecombs, int kmerSize, float windowPct) {
        boolean order = rng.nextBoolean();
        String[] seqs = order ? new String[] { seq1, seq2 } : new String[] { seq2, seq1 };
        int[] parents = order ? new int[] { 1, 2 } : new int[] { 2, 1 };

        Map<String, Integer> kmersA = kmerize(seqs[0], kmerSize);
        Map<String, Integer> kmersB = kmerize(seqs[1], kmerSize);

        Map<String, Pair<Integer, Integer>> kmersO = overlap(kmersA, kmersB, windowPct);
        List<String> kmers = new ArrayList<>(kmersO.keySet());

        if (kmers.size() <= numRecombs) { return null; }

        List<Pair<Integer, Integer>> recombPositions = new ArrayList<>();
        recombPositions.add(new Pair<>(0, 0));

        for (int r = 0; r < numRecombs; r++) {
            Pair<Integer, Integer> p;
            do {
                int index = rng.nextInt(kmers.size());

                p = kmersO.get(kmers.get(index));
            } while (p.getFirst() < recombPositions.get(recombPositions.size() - 1).getFirst() || p.getSecond() < recombPositions.get(recombPositions.size() - 1).getSecond());

            recombPositions.add(p);
        }

        recombPositions.add(new Pair<>(seqs[0].length() - 1, seqs[1].length() - 1));
        recombPositions.sort(Comparator.comparing(Pair::getFirst));

        List<String> pieces = new ArrayList<>();
        List<Integer> switches = new ArrayList<>();

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
            if (kmersA.containsKey(kmer) &&
                kmersB.containsKey(kmer) &&
                Math.abs(kmersA.get(kmer) - kmersB.get(kmer)) < (int) (windowPct*((float) kmersA.get(kmer)))) {
                kmersO.put(kmer, new Pair<>(kmersA.get(kmer), kmersB.get(kmer)));
            }
        }

        return kmersO;
    }
}
