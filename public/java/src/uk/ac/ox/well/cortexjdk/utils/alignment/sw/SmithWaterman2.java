package uk.ac.ox.well.cortexjdk.utils.alignment.sw;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * public domain snippet for SmithWaterman alignment
 * Smith, Temple F., and Michael S. Waterman.
 * "Identification of common molecular subsequences."
 * Journal of molecular biology 147.1 (1981): 195-197.
 *
 * author: yamule (https://github.com/yamule/smithwaterman)
 * usage ===
 *	SmithWaterman swold = new SmithWaterman();
 *	SWResult res = swold.align("EEEEMDQNNSLPPYAGGTWRYII","IIIIMDQNNSPPYAQGGTWRYEE");
 *	System.out.println("score: "+res.score);
 *	System.out.println(SmithWaterman.listToString(res.qseq));
 *	System.out.println(SmithWaterman.listToString(res.sseq));
 *
 * output ===
 * score: 73
 * ----EEEEMDQNNSLPPYA-GGTWRYII--
 * IIII----MDQNNS-PPYAQGGTWRY--EE
 *
 */
public class SmithWaterman2 {
    public static int TYPE_MATCH = 0;
    public static int TYPE_GAPINCOL = 1;
    public static int TYPE_GAPINROW = 2;
    SWCell[][] dpMat;
    //ScoringMatrix smat = new BLOSUM62();
    ScoringMatrix smat = new EDNAFULL();
    double penalO = 10;
    double penalE = 0.5;


    public SWResult align(String qq,String ss){
        String qq2 = qq.replaceAll("^[\\s]*>[^\\r\\n]*[\\r\\n]+","");
        String ss2 = ss.replaceAll("^[\\s]*>[^\\r\\n]*[\\r\\n]+","");

        String name1 = "seq1";
        String name2= "seq2";
        Pattern npat = Pattern.compile(">([^\\s]+)");
        Matcher mat = npat.matcher(qq);
        if(mat.find()){
            name1 = mat.group(1);
        }
        mat = npat.matcher(ss);
        if(mat.find()){
            name2 = mat.group(1);
        }

        Sequence q = new Sequence();
        Sequence s = new Sequence();
        q.name = name1;
        s.name = name2;
        q.add(qq2);
        s.add(ss2);
        return align(q,s);


    }

    public String[] getAlignment(String qq, String ss) {
        SWResult res = align(qq, ss);

        StringBuilder qqb = new StringBuilder();
        for (Character c : res.qseq) {
            qqb.append(c);
        }

        StringBuilder ssb = new StringBuilder();
        for (Character c : res.sseq) {
            ssb.append(c);
        }

        return new String[] { qqb.toString(), ssb.toString() };
    }

    public SWResult align(Sequence qq,Sequence ss){
        //Prepare letters for alignment.
        //Empty objects are removed and unknown letters are changed to 'X'.
        //Please see the filter method for exact process.
        ArrayList<Character> q = smat.filter(qq.seq);
        ArrayList<Character> s = smat.filter(ss.seq);
        int w = q.size()+1;
        int h = s.size()+1;

        dpMat = new SWCell[w][h];

        for(int xx = 0;xx < w;xx++){
            for(int yy = 0;yy < h;yy++){
                dpMat[xx][yy] = new SWCell();
            }
        }
        for(int xx = 0;xx < w;xx++){

            for(int ii = 0;ii < 3;ii++){
                dpMat[xx][0].setScoreAt(ii,0);
            }
        }
        for(int yy = 0;yy < h;yy++){
            for(int ii = 0;ii < 3;ii++){
                dpMat[0][yy].setScoreAt(ii,0);
            }
        }
        for(int xx = 1;xx <  w;xx++){
            for(int yy = 1;yy <  h;yy++){
                double fm = dpMat[xx-1][yy-1].getScoreAt(TYPE_MATCH) + smat.getScore(q.get(xx-1), s.get(yy-1));
                double fc = dpMat[xx-1][yy-1].getScoreAt(TYPE_GAPINCOL) + smat.getScore(q.get(xx-1), s.get(yy-1));
                double fr = dpMat[xx-1][yy-1].getScoreAt(TYPE_GAPINROW) + smat.getScore(q.get(xx-1), s.get(yy-1));
                SWCell currentcell = dpMat[xx][yy];
                fm = Math.max(0,fm);
                fc = Math.max(0,fc);
                fr = Math.max(0,fr);

                if(fm >= fr){
                    if(fm >= fc){
                        currentcell.setScoreAt(TYPE_MATCH,fm);
                        currentcell.setPrevTypeAt(TYPE_MATCH, TYPE_MATCH);
                    }else{
                        currentcell.setScoreAt(TYPE_MATCH,fc);
                        currentcell.setPrevTypeAt(TYPE_MATCH, TYPE_GAPINCOL);
                    }
                }else{
                    if(fr > fc){
                        currentcell.setScoreAt(TYPE_MATCH,fr);
                        currentcell.setPrevTypeAt(TYPE_MATCH, TYPE_GAPINROW);

                    }else{
                        currentcell.setScoreAt(TYPE_MATCH,fc);
                        currentcell.setPrevTypeAt(TYPE_MATCH, TYPE_GAPINCOL);
                    }
                }


                //Gap in col
                fm = dpMat[xx][yy-1].getScoreAt(TYPE_MATCH)-this.penalO;
                fc = dpMat[xx][yy-1].getScoreAt(TYPE_GAPINCOL)-this.penalE;
                fr = dpMat[xx][yy-1].getScoreAt(TYPE_GAPINROW)-this.penalO;// not used

                fm = Math.max(0,fm);
                fc = Math.max(0,fc);
                fr = Math.max(0,fr);


                if(fm >= fc){
                    currentcell.setScoreAt(TYPE_GAPINCOL, fm);
                    currentcell.setPrevTypeAt(TYPE_GAPINCOL, TYPE_MATCH);
                }else{
                    currentcell.setScoreAt(TYPE_GAPINCOL, fc);
                    currentcell.setPrevTypeAt(TYPE_GAPINCOL, TYPE_GAPINCOL);
                }

                //Gap in row
                fm = dpMat[xx-1][yy].getScoreAt(TYPE_MATCH)-this.penalO;
                fc = dpMat[xx-1][yy].getScoreAt(TYPE_GAPINCOL)-this.penalO;// not used
                fr = dpMat[xx-1][yy].getScoreAt(TYPE_GAPINROW)-this.penalE;

                fm = Math.max(0,fm);
                fc = Math.max(0,fc);
                fr = Math.max(0,fr);

                if(fm > fr){
                    currentcell.setScoreAt(TYPE_GAPINROW, fm);
                    currentcell.setPrevTypeAt(TYPE_GAPINROW, TYPE_MATCH);
                }else{
                    currentcell.setScoreAt(TYPE_GAPINROW, fr);
                    currentcell.setPrevTypeAt(TYPE_GAPINROW, TYPE_GAPINROW);
                }
            }
        }

        int maxpos[] = getMaxScorePos();
        int cx = maxpos[0];
        int cy = maxpos[1];
        int sx = cx;
        int sy = cy;
        int lx = cx;
        int ly = cy;
        if(cx == -1){//very short sequences which do not have positive value match between two.
            return new SWResult();
        }



        //backtrackpart
        double cs = dpMat[cx][cy].getScoreAt(TYPE_MATCH);
        double maxscore = cs;
        int ct = TYPE_MATCH;
        int cpt = dpMat[cx][cy].getPrevTypeAt(TYPE_MATCH);
        ArrayList<Character> qchar = new ArrayList<>();
        ArrayList<Character> schar = new ArrayList<>();
        qchar.add(q.get(cx-1));
        schar.add(s.get(cy-1));
        while(cs > 0){
            if(ct == TYPE_MATCH){
                cx--;
                cy--;
            }else if(ct == TYPE_GAPINCOL){
                cy--;
            }else if(ct == TYPE_GAPINROW){
                cx--;
            }
            ct = cpt;
            cs = dpMat[cx][cy].getScoreAt(ct);
            cpt = dpMat[cx][cy].getPrevTypeAt(ct);
            if(cs <= 0){
                break;
            }
            if(ct == TYPE_MATCH){
                qchar.add(q.get(cx-1));
                schar.add(s.get(cy-1));
                lx = cx;
                ly = cy;
            }else if(ct == TYPE_GAPINCOL){
                qchar.add('-');
                schar.add(s.get(cy-1));
                ly = cy;
            }else if(ct == TYPE_GAPINROW){
                qchar.add(q.get(cx-1));
                schar.add('-');
                lx = cx;
            }
        }


        //add unaligned nterminal part
        for(int xx = lx-1;xx > 0;xx--){
            qchar.add(q.get(xx-1));
            schar.add('-');
        }
        for(int yy = ly-1;yy > 0;yy--){
            schar.add(s.get(yy-1));
            qchar.add('-');
        }

        Collections.reverse(qchar);
        Collections.reverse(schar);


        //add unaligned cterminal part
        for(int xx = sx+1;xx < dpMat.length;xx++){
            qchar.add(q.get(xx-1));
            schar.add('-');
        }
        for(int yy = sy+1;yy < dpMat[0].length;yy++){
            schar.add(s.get(yy-1));
            qchar.add('-');
        }


        return new SWResult(qchar,schar,maxscore);

    }




    /**
     * Returns x y position of the cell which has maximum scores.
     * If there are multiple cells with the same score, the cell which has lowest x and lowest y is selected.
     * @return
     */
    public int[] getMaxScorePos(){
        int w = dpMat.length;
        int h = dpMat[0].length;
        int ret[] = new int[2];
        ret[0] = -1;
        ret[1] = -1;
        double maxscore = 0;
        for(int xx = 0;xx < w;xx++){
            for(int yy = 0;yy < h;yy++){
                for(int i = 0;i < 3;i++){
                    double sc = dpMat[xx][yy].getScoreAt(i);
                    if(maxscore < sc){
                        maxscore = sc;
                        ret[0] = xx;
                        ret[1] = yy;
                    }
                }
            }
        }
        return ret;
    }



    /**
     * Changes ArrayList<Character> to String.
     * @param al
     * @return
     */
    public static String listToString(ArrayList<Character> al){
        StringBuffer ret = new StringBuffer();
        for(Character c:al){
            ret.append(c);
        }
        return ret.toString();
    }

//    public static void main(String[] args){
//        for(int ii = 0;ii < args.length;ii++){
//            System.out.println(ii+" "+args[ii]);
//        }
//        ArrayList<Sequence> seq1 = Sequence.loadFasta(args[0]);
//        ArrayList<Sequence> seq2 = Sequence.loadFasta(args[1]);
//
//        SmithWaterman swold = new SmithWaterman2();
//        SWResult res = swold.align(seq1.get(0),seq2.get(0));
//        System.out.println("# score: "+res.score);
//        System.out.println(">"+seq1.get(0).name);
//        System.out.println(listToString(res.qseq));
//        System.out.println(">"+seq2.get(0).name);
//        System.out.println(listToString(res.sseq));
//
//		/*
//		// for debug
//		SmithWaterman swold = new SmithWaterman();
//		SWResult res = swold.align("MTPPPPGRAAPSAPRARVPGPPARLGLPLRLRLLLLLWAAAASAQGHLRSGPRIFAVWKG","PGPGNPSPMSLSPAWPGHPDQPLPREQMTSPAPPRIITSATADPEGTETALAGDTSDGLA");
//		System.out.println("score: "+res.score);
//		System.out.println(SmithWaterman.listToString(res.qseq));
//		System.out.println(SmithWaterman.listToString(res.sseq));
//		*/
//		/*
//		BLOSUM62 bl = new BLOSUM62();
//		SmithWaterman swold = new SmithWaterman();
//		ArrayList<Sequence> seqs = Sequence.loadFasta("C:\\Users\\kimidori\\Desktop\\TBM\\swold\\testfas.txt");
//		for(Sequence s:seqs){
//			System.out.println(">"+s.name+" "+s.desc);
//			for(Character c:s.seq){
//				System.out.print(String.valueOf(c));
//			}
//			System.out.println();
//		}
//		SWResult res = swold.align(seqs.get(0),seqs.get(1));
//		System.out.println(res.score);
//		System.out.println(seqs.get(0).name);
//		System.out.println(listToString(res.qseq));
//		System.out.println(seqs.get(1).name);
//		System.out.println(listToString(res.sseq));
//		*/
//    }
}



