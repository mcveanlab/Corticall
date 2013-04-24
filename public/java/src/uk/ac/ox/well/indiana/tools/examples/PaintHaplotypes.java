package uk.ac.ox.well.indiana.tools.examples;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationDiscrete;
import be.ac.ulg.montefiore.run.jahmm.OpdfDiscrete;
import be.ac.ulg.montefiore.run.jahmm.OpdfDiscreteFactory;
import be.ac.ulg.montefiore.run.jahmm.draw.GenericHmmDrawerDot;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class PaintHaplotypes extends Tool {
    @Argument(fullName="epsilon", shortName="e", doc="Error term")
    public Double EPSILON = 0.05;

    @Output
    public File out;

    private enum Copy {
        CORRECT, INCORRECT;

        public ObservationDiscrete<Copy> observation() {
            return new ObservationDiscrete<Copy>(this);
        }
    }

    private class StateMapping {
        public int donor;
        public int locus;

        public StateMapping(int d, int l) {
            donor = d;
            locus = l;
        }
    }

    private Map<Integer, StateMapping> getStateMap(String[] donors, String recipient) {
        Map<Integer, StateMapping> stateMap = new HashMap<Integer, StateMapping>();

        for (int i = 0, state = 0; i < donors.length; i++) {
            for (int j = 0; j < recipient.length(); j++, state++) {
                stateMap.put(state, new StateMapping(i, j));
            }
        }

        return stateMap;
    }

    private Hmm<ObservationDiscrete<Copy>> buildTrueHmm(String[] donors, String recipient, Map<Integer, StateMapping> stateMap) {
        double e = EPSILON / donors.length;

        Hmm<ObservationDiscrete<Copy>> hmm = new Hmm<ObservationDiscrete<Copy>>(donors.length * recipient.length(), new OpdfDiscreteFactory<Copy>(Copy.class));

        int state = 0;
        for (int i = 0; i < donors.length; i++) {
            for (int j = 0; j < recipient.length(); j++) {
                hmm.setPi(state, j == 0 ? 1.0 / donors.length : 0.0);

                double correctProb = donors[i].charAt(j) == recipient.charAt(j) ? 1.0 - EPSILON : EPSILON;

                hmm.setOpdf(state, new OpdfDiscrete<Copy>(Copy.class, new double[] { correctProb, 1.0 - correctProb }));

                state++;
            }
        }

        for (int state1 = 0; state1 < stateMap.size(); state1++) {
            for (int state2 = 0; state2 < stateMap.size(); state2++) {
                int donor1 = stateMap.get(state1).donor;
                int locus1 = stateMap.get(state1).locus;
                int donor2 = stateMap.get(state2).donor;
                int locus2 = stateMap.get(state2).locus;

                if (donor1 == donor2) { // same haplotype
                    if (locus2 == locus1 + 1) {
                        hmm.setAij(state1, state2, 1.0 - EPSILON);
                    } else {
                        hmm.setAij(state1, state2, 0.0);
                    }
                } else { // different haplotype
                    if (locus2 == locus1 + 1) {
                        hmm.setAij(state1, state2, EPSILON / (donors.length - 1));
                    } else {

                        hmm.setAij(state1, state2, 0.0);
                    }
                }
            }
        }

        return hmm;
    }

    private int[] relabelStatePath(int[] statePath, Map<Integer, StateMapping> stateMap) {
        int[] relabeledStatePath = new int[statePath.length];

        for (int i = 0; i < statePath.length; i++) {
            relabeledStatePath[i] = stateMap.get(statePath[i]).donor;
        }

        return relabeledStatePath;
    }

    @Override
    public void execute() {
        String[] donors = { "11100", "11000", "00111" };
        String[] recipients = { "00100", "11111", "01110" };

        List<ObservationDiscrete<Copy>> oseq = new ArrayList<ObservationDiscrete<Copy>>();
        for (int j = 0; j < recipients[0].length(); j++) {
            oseq.add(Copy.CORRECT.observation());
        }

        for (String donor : donors) {
            log.info("donor: {}", donor);
        }

        for (String recipient : recipients) {
            Map<Integer, StateMapping> stateMap = getStateMap(donors, recipient);

            Hmm<ObservationDiscrete<Copy>> hmm = buildTrueHmm(donors, recipient, stateMap);

            int[] statePath = hmm.mostLikelyStateSequence(oseq);
            int[] relabeledStatePath = relabelStatePath(statePath, stateMap);

            log.info("recip: {}, statePath: {}, relabeledStatePath: {}", recipient, statePath, relabeledStatePath);

            //log.info("hmm: {}", hmm);

            try {
                (new GenericHmmDrawerDot()).write(hmm, out.getAbsolutePath() + "." + recipient + ".dot");
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }


    }
}
