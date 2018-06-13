package uk.ac.ox.well.cortexjdk.utils.alignment.mosaic;

import htsjdk.samtools.util.StringUtil;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.utils.alignment.sw.SmithWaterman;

import java.util.*;

public class Tesserae {
    // constants
    private static final double SMALL = -1e32;

    // transition probabilities
    private double DEFAULT_DEL = 0.025;
    private double DEFAULT_EPS = 0.75;
    private double DEFAULT_REC = 0.0001;
    private double DEFAULT_TERM = 0.001;
    //private double DEFAULT_PMATCH = 0.0;

    // params
    private double del = DEFAULT_DEL;
    private double eps = DEFAULT_EPS;
    private double rho = DEFAULT_REC;
    private double term = DEFAULT_TERM;
    //private double pmatch = DEFAULT_PMATCH;

    private double ldel = Math.log(del);
    private double leps = Math.log(eps);
    private double lrho = Math.log(rho);
    private double lterm = Math.log(term);
    //private double lpmatch = Math.log(pmatch);

    private double piM = 0.75;
    private double piI = 1 - piM;
    private double mm = 1 - 2*del - rho - term;
    private double gm = 1 - eps - rho - term;
    private double dm = 1 - eps;

    private double lpiM = Math.log(piM);
    private double lpiI = Math.log(piI);
    private double lmm = Math.log(mm);
    private double lgm = Math.log(gm);
    private double ldm = Math.log(dm);

    private double[] emiss_gap_nt = { 0.2, 0.2, 0.2, 0.2, 0.2 };
    private double[][] emiss_match_nt = {
            {0.2, 0.2, 0.2, 0.2, 0.2},
            {0.2, 0.9, 0.05, 0.025, 0.025},
            {0.2, 0.05, 0.9, 0.025, 0.025},
            {0.2, 0.025, 0.025, 0.9, 0.05},
            {0.2, 0.025, 0.025, 0.05, 0.9},
    };

    // internal
    private int nseq = 0;
    private int maxl = 0;

    private double[][][] vt_m;
    private double[][][] vt_i;
    private double[][][] vt_d;
    private double[][][] tb_m;
    private double[][][] tb_i;
    private double[][][] tb_d;

    private int[] who_copy;
    private int[] maxpath_copy;
    private int[] maxpath_state;
    private int[] maxpath_pos;
    //private double[] ppsum_match;
    //private double[][] ppsum_state;

    private double tb_divisor;
    private double llk;
    private double combined_llk;

    private double[][] sm;
    private double[][] lsm;
    private double[] si;
    private double[] lsi;

    private List<Triple<String, String, Pair<Integer, Integer>>> path;
    private String editTrack;

    public Tesserae() { }

    public Tesserae(double del, double eps, double rho, double term) {
        this.del = del;
        this.eps = eps;
        this.rho = rho;
        this.term = term;
    }

    //public List<Triple<String, Pair<Integer, Integer>, String>> align(String query, Map<String, String> targets) {
    public List<Triple<String, String, Pair<Integer, Integer>>> align(String query, Map<String, String> targets) {
        initialize(query, targets);

        Map<String, String> panel = new LinkedHashMap<>();
        panel.put("query", query);
        panel.putAll(targets);

        return alignAll(panel);
    }

    private void initialize(String query, Map<String, String> targets) {
        // params
        ldel = Math.log(del);
        leps = Math.log(eps);
        lrho = Math.log(rho);
        lterm = Math.log(term);

        piM = 0.75;
        piI = 1 - piM;
        mm = 1 - 2*del - rho - term;
        gm = 1 - eps - rho - term;
        dm = 1 - eps;

        lpiM = Math.log(piM);
        lpiI = Math.log(piI);
        lmm = Math.log(mm);
        lgm = Math.log(gm);
        ldm = Math.log(dm);

        this.nseq = targets.size() + 1;
        this.maxl = getMaxLength(query, targets);

        vt_m = new double[nseq+1][maxl+1][maxl+1];
        vt_i = new double[nseq+1][maxl+1][maxl+1];
        vt_d = new double[nseq+1][maxl+1][maxl+1];
        tb_m = new double[nseq+1][maxl+1][maxl+1];
        tb_i = new double[nseq+1][maxl+1][maxl+1];
        tb_d = new double[nseq+1][maxl+1][maxl+1];

        who_copy = new int[nseq+1];
        Arrays.fill(who_copy, 1);
        who_copy[1] = 0;

        maxpath_copy = new int[2*maxl+1];
        maxpath_state = new int[2*maxl+1];
        maxpath_pos = new int[2*maxl+1];
        //ppsum_match = new double[nseq+1];
        //ppsum_state = new double[maxl+1][3];

        int tmp = (int) Math.log(maxl);
        tb_divisor = Math.pow(10.0, (double) tmp);

        combined_llk = 0.0;

        sm = new double[5][5];
        lsm = new double[5][5];

        for (int i = 0; i < sm.length; i++) {
            for (int j = 0; j < sm[0].length; j++) {
                sm[i][j] = emiss_match_nt[i][j];
                lsm[i][j] = Math.log(emiss_match_nt[i][j]);
            }
        }

        si = new double[5];
        lsi = new double[5];

        for (int i = 0; i < si.length; i++) {
            si[i] = emiss_gap_nt[i];
            lsi[i] = Math.log(emiss_gap_nt[i]);
        }

        path = new ArrayList<>();
    }

    private int getMaxLength(String query, Map<String, String> targets) {
        int maxLength = query.length();
        for (String target : targets.values()) {
            maxLength = Math.max(maxLength, target.length());
        }

        return maxLength;
    }

    private List<Triple<String, String, Pair<Integer, Integer>>> alignAll(Map<String, String> panel) {
        String query = panel.values().iterator().next();

        int l1 = query.length();
        double sizeL = 0.0;
        for (String target : panel.values()) {
            sizeL += target.length();
        }
        sizeL -= query.length();

        double lsizeL = Math.log(sizeL);

        int seq = 1;
        for (String target : panel.values()) {
            if (who_copy[seq] == 1) {
                int l2 = target.length();

                for (int pos_target = 0; pos_target <= l1; pos_target++) {
                    for (int pos_seq = 0; pos_seq <= l2; pos_seq++) {
                        vt_m[seq][pos_seq][pos_target] = SMALL;
                        vt_i[seq][pos_seq][pos_target] = SMALL;
                        vt_d[seq][pos_seq][pos_target] = SMALL;
                        tb_m[seq][pos_seq][pos_target] = 0;
                        tb_i[seq][pos_seq][pos_target] = 0;
                        tb_d[seq][pos_seq][pos_target] = 0;
                    }
                }
            }

            seq++;
        }

        int who_max = 0, state_max = 0, pos_max = 0, who_max_n = 0, state_max_n = 0, pos_max_n = 0;
        double max_r = SMALL;

        seq = 1;
        for (String target : panel.values()) {
            if (who_copy[seq] == 1) {
                for (int pos_seq = 1; pos_seq <= target.length(); pos_seq++) {
                    vt_m[seq][pos_seq][1] = lpiM - lsizeL;
                    vt_m[seq][pos_seq][1] += lsm[convert(query.charAt(0))][convert(target.charAt(pos_seq-1))];

                    vt_i[seq][pos_seq][1] = lpiI - lsizeL;
                    vt_i[seq][pos_seq][1] += lsi[convert(query.charAt(0))];

                    if (pos_seq > 0) {
                        vt_d[seq][pos_seq][1] = vt_m[seq][pos_seq - 1][1] + ldel;
                        tb_d[seq][pos_seq][1] = seq * 10 + 1 + (pos_seq - 1) / tb_divisor;

                        if ((vt_d[seq][pos_seq - 1][1] + leps) > vt_d[seq][pos_seq][1]) {
                            vt_d[seq][pos_seq][1] = vt_d[seq][pos_seq - 1][1] + leps;
                            tb_d[seq][pos_seq][1] = seq * 10 + 3 + (pos_seq - 1) / tb_divisor;
                        }
                    }

                    if (vt_m[seq][pos_seq][1] > max_r) {
                        max_r = vt_m[seq][pos_seq][1];
                        who_max = seq;
                        state_max = 1;
                        pos_max = pos_seq;
                    }
                    if (vt_i[seq][pos_seq][1] > max_r) {
                        max_r = vt_i[seq][pos_seq][1];
                        who_max = seq;
                        state_max = 2;
                        pos_max = pos_seq;
                    }
                }
            }

            seq++;
        }

        for (int pos_target = 2; pos_target <= l1; pos_target++) {
            double max_rn = SMALL + max_r;
            seq = 1;
            for (String target : panel.values()) {
                if (who_copy[seq] == 1) {
                    for (int pos_seq = 1; pos_seq <= target.length(); pos_seq++) {
                        // Match
                        vt_m[seq][pos_seq][pos_target] = max_r + lrho + lpiM - lsizeL;
                        tb_m[seq][pos_seq][pos_target] = who_max*10 + state_max + pos_max/tb_divisor;

                        // Compare to MM
                        if ((vt_m[seq][pos_seq-1][pos_target-1] + lmm) > vt_m[seq][pos_seq][pos_target]) {
                            vt_m[seq][pos_seq][pos_target] = vt_m[seq][pos_seq-1][pos_target-1] + lmm;
                            tb_m[seq][pos_seq][pos_target] = seq*10 + 1 + (pos_seq-1)/tb_divisor;
                        }

                        // Compare to IM
                        if ((vt_i[seq][pos_seq-1][pos_target-1] + lgm) > vt_m[seq][pos_seq][pos_target]) {
                            vt_m[seq][pos_seq][pos_target] = vt_i[seq][pos_seq-1][pos_target-1] + lgm;
                            tb_m[seq][pos_seq][pos_target] = seq*10 + 2 + (pos_seq-1)/tb_divisor;
                        }

                        // Compare to DM
                        if ((vt_d[seq][pos_seq-1][pos_target-1] + ldm) > vt_m[seq][pos_seq][pos_target]) {
                            vt_m[seq][pos_seq][pos_target] = vt_d[seq][pos_seq-1][pos_target-1] + ldm;
                            tb_m[seq][pos_seq][pos_target] = seq*10 + 3 + (pos_seq-1)/tb_divisor;
                        }

                        // Add in state match
                        vt_m[seq][pos_seq][pos_target] += lsm[convert(query.charAt(pos_target-1))][convert(target.charAt(pos_seq-1))];

                        // Insert
                        vt_i[seq][pos_seq][pos_target] = max_r + lrho + lpiI - lsizeL;
                        tb_i[seq][pos_seq][pos_target] = who_max*10 + state_max + pos_max/tb_divisor;

                        // Compare to MI
                        if ((vt_m[seq][pos_seq][pos_target-1] + ldel) > vt_i[seq][pos_seq][pos_target]) {
                            vt_i[seq][pos_seq][pos_target] = vt_m[seq][pos_seq][pos_target-1] + ldel;
                            tb_i[seq][pos_seq][pos_target] = seq*10 + 1 + pos_seq/tb_divisor;
                        }

                        // Compare to II
                        if ((vt_i[seq][pos_seq][pos_target-1] + leps) > vt_i[seq][pos_seq][pos_target]) {
                            vt_i[seq][pos_seq][pos_target] = vt_i[seq][pos_seq][pos_target-1] + leps;
                            tb_i[seq][pos_seq][pos_target] = seq*10 + 2 + pos_seq/tb_divisor;
                        }

                        // Add in state insert
                        vt_i[seq][pos_seq][pos_target] += lsi[convert(query.charAt(pos_target-1))];


                        // Delete
                        if (pos_target < l1 && pos_seq > 1) {
                            // Initialize with match
                            vt_d[seq][pos_seq][pos_target] = vt_m[seq][pos_seq-1][pos_target] + ldel;
                            tb_d[seq][pos_seq][pos_target] = seq*10 + 1 + (pos_seq-1)/tb_divisor;

                            // Compare to DD
                            if ((vt_d[seq][pos_seq-1][pos_target] + leps) > vt_d[seq][pos_seq][pos_target]) {
                                vt_d[seq][pos_seq][pos_target] = vt_d[seq][pos_seq-1][pos_target] + leps;
                                tb_d[seq][pos_seq][pos_target] = seq*10 + 3 + (pos_seq-1)/tb_divisor;
                            }
                        }

                        if (vt_m[seq][pos_seq][pos_target] > max_rn) {
                            max_rn = vt_m[seq][pos_seq][pos_target];
                            who_max_n = seq;
                            state_max_n = 1;
                            pos_max_n = pos_seq;
                        }
                        if (vt_i[seq][pos_seq][pos_target] > max_rn) {
                            max_rn = vt_i[seq][pos_seq][pos_target];
                            who_max_n = seq;
                            state_max_n = 2;
                            pos_max_n = pos_seq;
                        }
                    }
                }

                seq++;
            }

            max_r = max_rn;
            who_max = who_max_n;
            state_max = state_max_n;
            pos_max = pos_max_n;
        }

        llk = max_r + lterm;
        combined_llk += max_r + lterm;

        int cp = 2*maxl;
        maxpath_copy[cp] = who_max;
        maxpath_state[cp] = state_max;
        maxpath_pos[cp] = pos_max;

        int pos_target = l1;
        int who_next = 0, state_next = 0, pos_next = 0;
        while (pos_target >= 1) {
            if (state_max == 1) {
                who_next = (int) (tb_m[who_max][pos_max][pos_target]) / 10;
                state_next = (int) (tb_m[who_max][pos_max][pos_target] - who_next*10);
                pos_next = (int) ((tb_m[who_max][pos_max][pos_target] - who_next*10 - state_next)*tb_divisor + 1e-6);
            } else if (state_max == 2) {
                who_next = (int) (tb_i[who_max][pos_max][pos_target]) / 10;
                state_next = (int) (tb_i[who_max][pos_max][pos_target] - who_next*10);
                pos_next = (int) ((tb_i[who_max][pos_max][pos_target] - who_next*10-state_next)*tb_divisor + 1e-6);
            } else if (state_max == 3) {
                who_next = (int) (tb_d[who_max][pos_max][pos_target]) / 10;
                state_next = (int) (tb_d[who_max][pos_max][pos_target] - who_next*10);
                pos_next = (int) ((tb_d[who_max][pos_max][pos_target] - who_next*10-state_next)*tb_divisor + 1e-6);
            }

            cp--;

            maxpath_copy[cp] = who_next;
            maxpath_state[cp] = state_next;
            maxpath_pos[cp] = pos_next;

            who_max = who_next;
            state_max = state_next;
            pos_max = pos_next;

            if (maxpath_state[cp+1] != 3) {
                pos_target--;
            }
        }

        cp++;
        who_copy[1] = 1;

        // Prepare target sequence
        List<Pair<String, String>> seqs = new ArrayList<>();
        for (String seqName : panel.keySet()) {
            seqs.add(new Pair<>(seqName, panel.get(seqName)));
        }

        StringBuilder sb = new StringBuilder();
        int posStart = -1;
        int posEnd = -1;

        int i = 0;
        for (i = cp, pos_target = 1; i <= 2*maxl; i++) {
            if (maxpath_state[i] == 3) { sb.append("-"); }
            else {
                if (posStart == -1) {
                    posStart = pos_target-1;
                }
                posEnd = pos_target-1;

                sb.append(seqs.get(0).getValue().charAt(pos_target-1));
                pos_target++;
            }
        }

        path.add(Triple.of(seqs.get(0).getFirst(), sb.toString(), Pair.create(posStart, posEnd)));

        // Prepare matching track
        sb = new StringBuilder();
        for (i = cp, pos_target = 1; i <= 2*maxl; i++) {
            if (maxpath_state[i] == 1) {
                if (seqs.get(0).getSecond().charAt(pos_target - 1) == seqs.get(maxpath_copy[i]-1).getSecond().charAt(maxpath_pos[i]-1)) {
                    sb.append("|");
                } else {
                    sb.append(" ");
                }

                pos_target++;
            } else if (maxpath_state[i] == 2) {
                pos_target++;
                sb.append("^");
            } else {
                sb.append("~");
            }
        }

        editTrack = sb.toString();

        // Prepare copying tracks
        String currentTrack = seqs.get(maxpath_copy[cp]-1).getFirst();
        sb = new StringBuilder();
        posStart = -1;
        posEnd = -1;
        int lastKnownPos = -1;

        boolean uppercase = true;
        for (i = cp; i <= 2*maxl; i++) {
            if (i > cp && maxpath_copy[i] == maxpath_copy[i-1] && Math.abs(maxpath_pos[i] - maxpath_pos[i - 1]) > 1 || maxpath_pos[i] == lastKnownPos + 1) {
                path.add(Triple.of(currentTrack, sb.toString(), Pair.create(posStart, posEnd)));
                uppercase = !uppercase;
                lastKnownPos = maxpath_pos[i - 1];

                if (posStart != posEnd) {
                    posStart = maxpath_pos[i] - 1;
                    posEnd = maxpath_pos[i] - 1;
                }

                currentTrack = seqs.get(maxpath_copy[i]-1).getFirst();
                sb = new StringBuilder();
                sb.append(StringUtil.repeatCharNTimes(' ', i - cp));
            }

            if (i > cp && maxpath_copy[i] != maxpath_copy[i-1]) {
                path.add(Triple.of(currentTrack, sb.toString(), Pair.create(posStart, posEnd)));
                uppercase = true;

                if (posStart != posEnd) {
                    posStart = maxpath_pos[i] - 1;
                    posEnd = maxpath_pos[i] - 1;
                }

                currentTrack = seqs.get(maxpath_copy[i]-1).getFirst();
                sb = new StringBuilder();
                sb.append(StringUtil.repeatCharNTimes(' ', i - cp));
            }

            if (maxpath_state[i] == 2) {
                sb.append("-");
            } else {
                char c = seqs.get(maxpath_copy[i] - 1).getSecond().charAt(maxpath_pos[i] - 1);
                c = uppercase ? Character.toUpperCase(c) : Character.toLowerCase(c);

                if (posStart == -1) {
                    posStart = maxpath_pos[i] - 1;
                }
                posEnd = maxpath_pos[i] - 1;

                sb.append(c);
            }
        }

        path.add(Triple.of(currentTrack, sb.toString(), Pair.create(posStart, posEnd)));

        return path;
    }

    private int convert(char c) {
        switch (c) {
            case 'A': return 3;
            case 'C': return 2;
            case 'G': return 4;
            case 'T': return 1;
        }

        return 0;
    }

    public double getMaximumLogLikelihood() {
        return llk;
    }

    @Override
    public String toString() {
        int maxNameLength = 0;
        for (int i = 0; i < path.size(); i++) {
            String name = String.format("%s (%d-%d)", path.get(i).getLeft(), path.get(i).getRight().getFirst(), path.get(i).getRight().getSecond());
            maxNameLength = Math.max(maxNameLength, name.length());
        }
        String format = "%" + maxNameLength + "s";

        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < path.size(); i++) {
            String name = String.format("%s (%d-%d)", path.get(i).getLeft(), path.get(i).getRight().getFirst(), path.get(i).getRight().getSecond());

            sb.append(String.format(format, name))
              .append(" ")
              .append(path.get(i).getMiddle())
              .append("\n");

            if (i == 0) {
                sb.append(String.format(format, " "))
                  .append(" ")
                  .append(editTrack)
                  .append("\n");
            }
        }

        sb.append("\n")
          .append("Mllk: ")
          .append(llk)
          .append("\n");

        return sb.toString();
    }
}
