package uk.ac.ox.well.cortexjdk.utils.alignment.sw;

import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;

import java.util.ArrayList;
import java.util.HashSet;

class EDNAFULL implements ScoringMatrix{
    String lines[] = {
            //ftp://ftp.ncbi.nih.gov/blast/matrices/NUC.4.4
            "    A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N   *",
            "A   5  -4  -4  -4  -4   1   1  -4  -4   1  -4  -1  -1  -1  -2  -4",
            "T  -4   5  -4  -4  -4   1  -4   1   1  -4  -1  -4  -1  -1  -2  -4",
            "G  -4  -4   5  -4   1  -4   1  -4   1  -4  -1  -1  -4  -1  -2  -4",
            "C  -4  -4  -4   5   1  -4  -4   1  -4   1  -1  -1  -1  -4  -2  -4",
            "S  -4  -4   1   1  -1  -4  -2  -2  -2  -2  -1  -1  -3  -3  -1  -4",
            "W   1   1  -4  -4  -4  -1  -2  -2  -2  -2  -3  -3  -1  -1  -1  -4",
            "R   1  -4   1  -4  -2  -2  -1  -4  -2  -2  -3  -1  -3  -1  -1  -4",
            "Y  -4   1  -4   1  -2  -2  -4  -1  -2  -2  -1  -3  -1  -3  -1  -4",
            "K  -4   1   1  -4  -2  -2  -2  -2  -1  -4  -1  -3  -3  -1  -1  -4",
            "M   1  -4  -4   1  -2  -2  -2  -2  -4  -1  -3  -1  -1  -3  -1  -4",
            "B  -4  -1  -1  -1  -1  -3  -3  -1  -1  -3  -1  -2  -2  -2  -1  -4",
            "V  -1  -4  -1  -1  -1  -3  -1  -3  -3  -1  -2  -1  -2  -2  -1  -4",
            "H  -1  -1  -4  -1  -3  -1  -3  -1  -3  -1  -2  -2  -1  -2  -1  -4",
            "D  -1  -1  -1  -4  -3  -1  -1  -3  -1  -3  -2  -2  -2  -1  -1  -4",
            "N  -2  -2  -2  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -4",
            "*  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4   1",
    };


    HashSet<Character> acceptable = new HashSet<>();
    int[][] scoreMat = new int[128][128];

    EDNAFULL(){
        for(int ii = 0;ii < 128;ii++){
            for(int jj = 0;jj < 128;jj++){
                scoreMat[ii][jj] = 10000;
            }
        }
        ArrayList<String> head = splitWithSpace(lines[0]);
        for(int ii = 1;ii < lines.length;ii++){
            ArrayList<String> pt = splitWithSpace(lines[ii]);
            String s = pt.get(0);
            char sc = s.charAt(0);
            acceptable.add(sc);
            for(int jj = 1;jj < pt.size();jj++){
                String q = head.get(jj-1);
                char qc = q.charAt(0);
                acceptable.add(qc);
                if(scoreMat[sc][qc] < 100){
                    throw new CortexJDKException(sc+"-"+qc+" duplicate pair?");
                }
                scoreMat[sc][qc] = Integer.parseInt(pt.get(jj));
            }
        }
    }


    /**
     * Change letters which are not compatible with this matrix.
     * @param al
     */
    public ArrayList<Character> filter(ArrayList<Character> al){
        ArrayList<Character> ret = new ArrayList<>();
        for(int ii = 0;ii < al.size();ii++){
            Character a = al.get(ii);
            if(a == null || a == Character.MIN_VALUE){
            }else{
                if(a.equals('-') || a.equals('.')){
                }else{
                    if(acceptable.contains(a)){
                        ret.add(a);
                    }else{
                        ret.add('X');
                    }
                }
            }
        }
        return ret;
    }



    //For discrepancy between open java & oracle java
    public ArrayList<String> splitWithSpace(String str){
        ArrayList<String> ret = new ArrayList<>();
        String head[] = str.replaceAll("^[\\s]+","").replaceAll("[\\s]+$","").split("[\\s]+");
        for(String s:head){
            if(s.length() > 0){
                ret.add(s);
            }
        }
        return ret;
    }

    public int getScore(char q,char s){
        return scoreMat[q][s];
    }
}
