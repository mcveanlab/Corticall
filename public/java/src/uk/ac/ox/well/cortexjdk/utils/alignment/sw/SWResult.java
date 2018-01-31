package uk.ac.ox.well.cortexjdk.utils.alignment.sw;

import java.util.ArrayList;

/**
 * Stores result of SmithWaterman Alignment.
 * @author kimidori
 */
public class SWResult {
    ArrayList<Character> qseq = new ArrayList<>();
    ArrayList<Character> sseq = new ArrayList<>();
    double score = -1;
    SWResult(){
    }
    SWResult(ArrayList<Character> qq,ArrayList<Character> ss,double s){
        qseq.addAll(qq);
        sseq.addAll(ss);
        score = s;
    }
}
