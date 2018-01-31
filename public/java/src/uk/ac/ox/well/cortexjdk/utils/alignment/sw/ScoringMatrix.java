package uk.ac.ox.well.cortexjdk.utils.alignment.sw;

import java.util.ArrayList;

interface ScoringMatrix{
    public int getScore(char q, char s);
    public ArrayList<Character> filter(ArrayList<Character> al);
}
