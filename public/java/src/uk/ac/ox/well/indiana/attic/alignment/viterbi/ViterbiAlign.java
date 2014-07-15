package uk.ac.ox.well.indiana.attic.alignment.viterbi;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleFactory3D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.DoubleMatrix3D;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ViterbiAlign extends Module {
    @Argument(fullName="seq", shortName="seq", doc="Name of file containing sequences in FASTA format")
    public FastaSequenceFile SEQ;

    @Argument(fullName="params", shortName="params", doc="Name of file containing input parameters")
    public File PARAMS;

    @Argument(fullName="target", shortName="target", doc="Group to analyze as target sequences")
    public ArrayList<String> TARGETS;

    @Argument(fullName="rec", shortName="rec", doc="Recombination probability")
    public Float REC = 0.0001f;

    @Argument(fullName="del", shortName="del", doc="Rate of moving to delete state")
    public Float DEL = 0.025f;

    @Argument(fullName="eps", shortName="eps", doc="Rate of extending in delete/insert state")
    public Float EPS = 0.750f;

    @Argument(fullName="term", shortName="term", doc="Rate of moving to termination state")
    public Float TERM = 0.001f;

    @Argument(fullName="pmatch", shortName="pmatch", doc="Probability of match")
    public Float PMATCH = 0.0f;

    @Argument(fullName="tag", shortName="tag", doc="Tag to attach to output files")
    public String TAG = "";

    @Argument(fullName="psum", shortName="psum", doc="Print out summed posteriors for match states in each sequence")
    public Boolean PSUM = false;

    @Argument(fullName="estimate", shortName="estimate", doc="Estimate parameters by EM with Pr{rec}=0")
    public Boolean ESTIMATE = false;

    private double SMALL = -1e32;
    private int NSTATE_AA = 25;
    private DoubleMatrix2D lsm = DoubleFactory2D.dense.make(NSTATE_AA, NSTATE_AA, 0.0);

    private void initialize() {

    }

    private class Matrices {
        private DoubleMatrix3D vt_m;
        private DoubleMatrix3D vt_i;
        private DoubleMatrix3D vt_d;

        private DoubleMatrix3D tb_m;
        private DoubleMatrix3D tb_i;
        private DoubleMatrix3D tb_d;

        private int maxpathCopy;
        private int maxpathState;
        private int maxpathPos;
        private double ppsumMatch;
        private double ppsumState;

        private double llk;

        private double tbDivisor;

        private double[] expectedTransitions;
        private double[] expectedEmissions;

        public void reset() {
            for (int s = 0; s < vt_d.slices(); s++) {
                for (int r = 0; r < vt_m.rows(); r++) {
                    for (int c = 0; c < vt_m.columns(); c++) {
                        vt_m.set(s, r, c, SMALL);
                        tb_m.set(s, r, c, 0);
                    }
                }
            }
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();

            sb.append("\n");
            sb.append("vt_m: " + vt_m.slices() + "x" + vt_m.rows() + "x" + vt_m.columns() + "\n");
            sb.append("vt_i: " + vt_i.slices() + "x" + vt_i.rows() + "x" + vt_i.columns() + "\n");
            sb.append("vt_d: " + vt_d.slices() + "x" + vt_d.rows() + "x" + vt_d.columns() + "\n");
            sb.append("tb_m: " + tb_m.slices() + "x" + tb_m.rows() + "x" + tb_m.columns() + "\n");
            sb.append("tb_i: " + tb_i.slices() + "x" + tb_i.rows() + "x" + tb_i.columns() + "\n");
            sb.append("tb_d: " + tb_d.slices() + "x" + tb_d.rows() + "x" + tb_d.columns() + "\n");

            return sb.toString();
        }
    }

    private Map<String, List<ReferenceSequence>> loadSequences() {
        Map<String, List<ReferenceSequence>> seqs = new HashMap<String, List<ReferenceSequence>>();
        seqs.put("target", new ArrayList<ReferenceSequence>());
        seqs.put("ref", new ArrayList<ReferenceSequence>());

        ReferenceSequence rseq;
        while ((rseq = SEQ.nextSequence()) != null) {
            String name = rseq.getName();

            boolean isTarget = false;
            for (String target : TARGETS) {
                if (name.contains(target)) {
                    seqs.get("target").add(rseq);

                    isTarget = true;
                    break;
                }
            }

            if (!isTarget) {
                seqs.get("ref").add(rseq);
            }
        }

        return seqs;
    }

    private Matrices allocateMatrices(Map<String, List<ReferenceSequence>> myData) {
        int maxl = 0;

        for (String group : myData.keySet()) {
            for (ReferenceSequence rseq : myData.get(group)) {
                if (rseq.length() > maxl) {
                    maxl = rseq.length();
                }
            }
        }

        int nseq = myData.get("ref").size();

        Matrices matrices = new Matrices();
        matrices.vt_m = DoubleFactory3D.dense.make(nseq, maxl, maxl);
        matrices.vt_i = DoubleFactory3D.dense.make(nseq, maxl, maxl);
        matrices.vt_d = DoubleFactory3D.dense.make(nseq, maxl, maxl);
        matrices.tb_m = DoubleFactory3D.dense.make(nseq, maxl, maxl);
        matrices.tb_i = DoubleFactory3D.dense.make(nseq, maxl, maxl);
        matrices.tb_d = DoubleFactory3D.dense.make(nseq, maxl, maxl);

        return matrices;
    }

    /*
        Ala     A 1      Alanine
        Arg     R 2      Arginine
        Asn     N 3      Asparagine
        Asp     D 4      Aspartic acid (Aspartate)
        Cys     C 5      Cysteine
        Gln     Q 6      Glutamine
        Glu     E 7      Glutamic acid (Glutamate)
        Gly     G 8      Glycine
        His     H 9      Histidine
        Ile     I 10     Isoleucine
        Leu     L 11     Leucine
        Lys     K 12     Lysine
        Met     M 13     Methionine
        Phe     F 14     Phenylalanine
        Pro     P 15     Proline
        Ser     S 16     Serine
        Thr     T 17     Threonine
        Trp     W 18     Tryptophan
        Tyr     Y 19     Tyrosine
        Val     V 20     Valine
        Asx     B 21     Aspartic acid or Asparagine
        Glx     Z 22     Glutamine or Glutamic acid.
        Xaa     X 23     Any amino acid.
        TERM      24     termination codon
     */

    private int[] recodeSequence(ReferenceSequence rseq) {
        String seq = new String(rseq.getBases());
        int[] recodedSeq = new int[rseq.length()];

        for (int i = 0; i < recodedSeq.length; i++) {
            switch (seq.charAt(i)) {
                case 'A': case 'a': recodedSeq[i] = 1; break; // Alanine
                case 'R': case 'r': recodedSeq[i] = 2; break; // Arginine
                case 'N': case 'n': recodedSeq[i] = 3; break; // Asparagine
                case 'D': case 'd': recodedSeq[i] = 4; break; // Aspartic acid (Aspartate)
                case 'C': case 'c': recodedSeq[i] = 5; break; // Cysteine
                case 'Q': case 'q': recodedSeq[i] = 6; break; // Glutamine
                case 'E': case 'e': recodedSeq[i] = 7; break; // Glutamic acid (Glutamate)
                case 'G': case 'g': recodedSeq[i] = 8; break; // Glycine
                case 'H': case 'h': recodedSeq[i] = 9; break; // Histidine
                case 'I': case 'i': recodedSeq[i] = 10; break; // Isoleucine
                case 'L': case 'l': recodedSeq[i] = 11; break; // Leucine
                case 'K': case 'k': recodedSeq[i] = 12; break; // Lysine
                case 'M': case 'm': recodedSeq[i] = 13; break; // Methionine
                case 'F': case 'f': recodedSeq[i] = 14; break; // Phenylalanine
                case 'P': case 'p': recodedSeq[i] = 15; break; // Proline
                case 'S': case 's': recodedSeq[i] = 16; break; // Serine
                case 'T': case 't': recodedSeq[i] = 17; break; // Threonine
                case 'W': case 'w': recodedSeq[i] = 18; break; // Tryptophan
                case 'Y': case 'y': recodedSeq[i] = 19; break; // Tyrosine
                case 'V': case 'v': recodedSeq[i] = 20; break; // Valine
                case 'B': case 'b': recodedSeq[i] = 21; break; // Aspartic acid or Asparagine
                case 'Z': case 'z': recodedSeq[i] = 22; break; // Glutamine or Glutamic acid.
                case 'X': case 'x': recodedSeq[i] = 23; break; // Any amino acid.
                case '.': case '*': recodedSeq[i] = 24; break; // termination codon
                default:
                    throw new IndianaException("Could not convert '" + seq.charAt(i) + "' to an amino acid code.");
            }
        }

        return recodedSeq;
    }

    private void kalignVt(Map<String, List<ReferenceSequence>> myData, Matrices myMatrices, ReferenceSequence target) {
        int sizeL = 0;
        for (ReferenceSequence rseq : myData.get("ref")) {
            sizeL += rseq.length();
        }

        double lsizeL = Math.log((double) sizeL);

        log.info("sizeL: {}", sizeL);
        log.info("lsizeL: {}", lsizeL);

        int l1 = target.length(); // length of target sequence

        myMatrices.reset();

        double piM = 0.75;
        double lpiM = Math.log(piM);

        // Initialize Viterbi matrices
        List<ReferenceSequence> refs = myData.get("ref");
        for (int seq = 0; seq < refs.size(); seq++) {
            //for (int pos_seq = 1; pos_seq <= l2; )
            /*Viterbi Ms*/
//            vt_m[seq][pos_seq][1] = (double) my_pars->lpiM-my_pars->lsizeL;
//            vt_m[seq][pos_seq][1] += my_pars->lsm[s1[1]][s2[pos_seq]];
//            vt_i[seq][pos_seq][1] = my_pars->lpiI-my_pars->lsizeL;
//            vt_i[seq][pos_seq][1] += my_pars->lsi[s1[1]];
            //myMatrices.vt_m.set(seq, );

            ReferenceSequence ref = refs.get(seq);
            int l2 = ref.length();
            int[] s2 = recodeSequence(ref);

            for (int pos_seq = 1; pos_seq <= l2; pos_seq++) {
                myMatrices.vt_m.set(seq, pos_seq, 1, lpiM - lsizeL);
                //myMatrices.vt_m.set(seq, pos_seq, 1, myMatrices.vt_m.get(seq, pos_seq, 1) +
            }
        }
    }

    @Override
    public void execute() {
        Map<String, List<ReferenceSequence>> myData = loadSequences();
        Matrices myMatrices = allocateMatrices(myData);

        log.info("matrices: {}", myMatrices);

        for (ReferenceSequence rseq : myData.get("target")) {
            kalignVt(myData, myMatrices, rseq);
        }
    }
}
