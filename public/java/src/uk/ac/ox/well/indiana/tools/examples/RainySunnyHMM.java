package uk.ac.ox.well.indiana.tools.examples;

import be.ac.ulg.montefiore.run.jahmm.*;
import uk.ac.ox.well.indiana.tools.Tool;

import java.util.*;

public class RainySunnyHMM extends Tool {
    private enum Task {
        WALK, SHOP, CLEAN;

        public ObservationDiscrete<Task> observation() {
            return new ObservationDiscrete<Task>(this);
        }
    }

    /**
     * Build two-state HMM (rainy, sunny)
     *
     * @return  HMM
     */
    private static Hmm<ObservationDiscrete<Task>> buildTrueHmm() {
        Hmm<ObservationDiscrete<Task>> hmm = new Hmm<ObservationDiscrete<Task>>(2, new OpdfDiscreteFactory<Task>(Task.class));

        hmm.setPi(0, 0.6);
        hmm.setPi(1, 0.4);

        hmm.setOpdf(0, new OpdfDiscrete<Task>(Task.class, new double[] { 0.1, 0.4, 0.5 }));
        hmm.setOpdf(1, new OpdfDiscrete<Task>(Task.class, new double[] { 0.6, 0.3, 0.1 }));

        hmm.setAij(0, 0, 0.7);
        hmm.setAij(0, 1, 0.3);
        hmm.setAij(1, 0, 0.4);
        hmm.setAij(1, 1, 0.6);

        return hmm;
    }

    @Override
    public void execute() {
        Hmm<ObservationDiscrete<Task>> hmm = buildTrueHmm();

        List<ObservationDiscrete<Task>> seq = new ArrayList<ObservationDiscrete<Task>>();
        seq.add(Task.WALK.observation());
        seq.add(Task.CLEAN.observation());
        seq.add(Task.SHOP.observation());

        ForwardBackwardCalculator fbc = new ForwardBackwardCalculator(seq, hmm, EnumSet.of(ForwardBackwardCalculator.Computation.ALPHA));

        log.info("Seq: {}, prob: {}", seq, fbc.probability());
    }
}
