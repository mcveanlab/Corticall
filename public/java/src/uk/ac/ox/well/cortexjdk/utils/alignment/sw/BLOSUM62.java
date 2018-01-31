package uk.ac.ox.well.cortexjdk.utils.alignment.sw;

import java.util.ArrayList;
import java.util.HashSet;

class BLOSUM62 implements ScoringMatrix{
    String lines[] = {
            //https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
            "   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *",
            "A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 ",
            "R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 ",
            "N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 ",
            "D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 ",
            "C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 ",
            "Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 ",
            "E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 ",
            "G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 ",
            "H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 ",
            "I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 ",
            "L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 ",
            "K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 ",
            "M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 ",
            "F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 ",
            "P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 ",
            "S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 ",
            "T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 ",
            "W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 ",
            "Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 ",
            "V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 ",
            "B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 ",
            "Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 ",
            "X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 ",
            "* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 "
    };


    HashSet<Character> acceptable = new HashSet<>();
    int[][] scoreMat = new int[128][128];

    BLOSUM62(){
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
                    System.out.println(sc+"-"+qc+" duplicate pair?");
                }
                scoreMat[sc][qc] = Integer.parseInt(pt.get(jj));
            }
        }
		/*
		ArrayList<Character> aa = new ArrayList<>(acceptable);
		for(Character k:aa){
			for(Character q:aa){
				int sc3 = scoreMat[q][k];
				int sc4 = scoreMat[k][q];
				if(sc3 == sc4){
					System.out.println(k+"-"+q+" OK "+sc3);
				}else{
					System.out.println(k+q+" is different, "+sc3+","+sc4+",");
				}
			}
		}
		*/

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
