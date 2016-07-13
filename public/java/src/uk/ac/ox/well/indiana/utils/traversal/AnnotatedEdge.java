package uk.ac.ox.well.indiana.utils.traversal;

import org.jgrapht.graph.DefaultWeightedEdge;

/**
 * Created by kiran on 18/06/2015.
 */
public class AnnotatedEdge extends DefaultWeightedEdge {
    private boolean[] isInColor = new boolean[10];

    public AnnotatedEdge() {}

    public AnnotatedEdge(boolean... cs) {
        for (int c = 0; c < cs.length; c++) {
            set(c, cs[c]);
        }
    }

    public AnnotatedEdge(int presentInColor) {
        setPresence(presentInColor);
    }

    public void setPresence(int c) {
        set(c, true);
    }

    public void setAbsence(int c) {
        set(c, false);
    }

    public void set(int c, boolean presence) {
        isInColor[c] = presence;
    }

    public boolean isPresent(int c) {
        return isInColor[c];
    }

    public boolean isAbsent(int c) { return !isInColor[c]; }

    public boolean[] getPresence() { return isInColor; }
}
