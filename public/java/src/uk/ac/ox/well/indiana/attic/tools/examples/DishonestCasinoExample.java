package uk.ac.ox.well.indiana.attic.tools.examples;

import be.ac.ulg.montefiore.run.jahmm.*;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchLearner;
import be.ac.ulg.montefiore.run.jahmm.toolbox.KullbackLeiblerDistanceCalculator;
import be.ac.ulg.montefiore.run.jahmm.toolbox.MarkovGenerator;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;

import java.util.ArrayList;
import java.util.List;

public class DishonestCasinoExample extends Module {
    @Argument(fullName="iterations", shortName="i", doc="Number of iterations for HMM refinement")
    public Integer ITERATIONS = 10;

    public enum Die {
        ONE, TWO, THREE, FOUR, FIVE, SIX;

        public ObservationDiscrete<Die> observation() {
            return new ObservationDiscrete<Die>(this);
        }
    }

    private static Hmm<ObservationDiscrete<Die>> buildTrueHmm() {
        Hmm<ObservationDiscrete<Die>> hmm = new Hmm<ObservationDiscrete<Die>>(2, new OpdfDiscreteFactory<Die>(Die.class));

        hmm.setPi(0, 1.0);
        hmm.setPi(1, 0.0);

        hmm.setOpdf(0, new OpdfDiscrete<Die>(Die.class, new double[] { 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0 }));
        hmm.setOpdf(1, new OpdfDiscrete<Die>(Die.class, new double[] { 1.0/10.0, 1.0/10.0, 1.0/10.0, 1.0/10.0, 1.0/10.0, 1.0/2.0 }));

        hmm.setAij(0, 0, 0.95);
        hmm.setAij(0, 1, 0.05);
        hmm.setAij(1, 0, 0.10);
        hmm.setAij(1, 1, 0.90);

        return hmm;
    }

    private static Hmm<ObservationDiscrete<Die>> buildInitialGuessHmm() {
        Hmm<ObservationDiscrete<Die>> hmm = new Hmm<ObservationDiscrete<Die>>(2, new OpdfDiscreteFactory<Die>(Die.class));

        hmm.setPi(0, 0.5);
        hmm.setPi(1, 0.5);

        hmm.setOpdf(0, new OpdfDiscrete<Die>(Die.class, new double[] { 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0 }));
        hmm.setOpdf(1, new OpdfDiscrete<Die>(Die.class, new double[] { 0.133, 0.133, 0.133, 0.133, 0.133, 1.0/3.0 }));

        hmm.setAij(0, 0, 0.90);
        hmm.setAij(0, 1, 0.10);
        hmm.setAij(1, 0, 0.10);
        hmm.setAij(1, 1, 0.90);

        return hmm;
    }

    private static <O extends Observation> List<List<O>> generateSequences(Hmm<O> hmm) {
        MarkovGenerator<O> mg = new MarkovGenerator<O>(hmm);

        List<List<O>> sequences = new ArrayList<List<O>>();

        for (int i = 0; i < 200; i++) {
            sequences.add(mg.observationSequence(100));
        }

        return sequences;
    }

    @Override
    public void execute() {
        // Initialize a hidden markov model with two hidden states (fair vs. loaded)
        Hmm<ObservationDiscrete<Die>> hmm = buildTrueHmm();

        // Generate a sequence of die rolls, switching between dies as specified by the HMM
        List<List<ObservationDiscrete<Die>>> sequences = generateSequences(hmm);

        // Now guess an initial HMM.  We assume there are two hidden states.
        Hmm<ObservationDiscrete<Die>> learnedHmm = buildInitialGuessHmm();

        BaumWelchLearner bwl = new BaumWelchLearner();
        KullbackLeiblerDistanceCalculator klc = new KullbackLeiblerDistanceCalculator();

        // Iterate and learn HMM parameters
        for (int i = 0; i < ITERATIONS; i++) {
            log.info("iteration {}: hmm distance: {}", i, klc.distance(hmm, learnedHmm));

            learnedHmm = bwl.iterate(learnedHmm, sequences);
        }

        log.info("initial hmm: {}", hmm);
        log.info("  final hmm: {}", learnedHmm);

        MarkovGenerator<ObservationDiscrete<Die>> mg = new MarkovGenerator<ObservationDiscrete<Die>>(hmm);
        List<ObservationDiscrete<Die>> obsSequence = new ArrayList<ObservationDiscrete<Die>>();
        List<Integer> stateSequence = new ArrayList<Integer>();

        for (int i = 0; i < 100; i++) {
            stateSequence.add(mg.stateNb());
            obsSequence.add(mg.observation());
        }

        ViterbiCalculator vc = new ViterbiCalculator(obsSequence, hmm);

        String seq = "";
        String states = "";
        String trueStates = "";
        for (int b = 0; b < sequences.get(0).size(); b++) {
            seq += sequences.get(0).get(b).value.ordinal() + 1;
            states += vc.stateSequence()[b];
            trueStates += stateSequence.get(b);
        }

        log.info(seq);
        log.info(states);
        log.info(trueStates);
    }
}
