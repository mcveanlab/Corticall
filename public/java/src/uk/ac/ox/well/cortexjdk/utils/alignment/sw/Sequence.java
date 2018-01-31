package uk.ac.ox.well.cortexjdk.utils.alignment.sw;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

class Sequence{
    String name;
    String desc;
    ArrayList<Character> seq = new ArrayList<Character>();
    public static ArrayList<Sequence> loadFasta(String filename){
        ArrayList<Sequence> ret = new ArrayList<>();
        try{
            BufferedReader br = new  BufferedReader(new FileReader(new File(filename)));
            String ln = null;
            Sequence currentseq = new Sequence();
            ret.add(currentseq);
            Pattern pat1 = Pattern.compile("^[\\s]*>[\\s]*([^\\s]*)");
            Pattern pat2 = Pattern.compile("^[\\s]*>[\\s]*([^\\s]+)[\\s]+([^\\r\\n]+)");
            while((ln = br.readLine()) != null){
                Matcher mat = pat1.matcher(ln);
                if(mat.find()){
                    String n = mat.group(1);
                    String d = "";
                    Matcher mat2 = pat2.matcher(ln);
                    if(mat2.find()){
                        d = mat2.group(2);
                    }
                    if(ret.size() == 1 && currentseq.seq.size() == 0){
                    }else{
                        currentseq = new Sequence();
                        ret.add(currentseq);
                    }
                    currentseq.name = n;
                    currentseq.desc = d;
                }else{
                    currentseq.add(ln);
                }
            }
        }catch(Exception exx){
            exx.printStackTrace();
        }
        return ret;
    }
    public void add(String s){
        String[] pt = s.replaceAll("[\\r\\n]","").split("");
        for(String pp:pt){
            if(pp.length() == 1){
                seq.add(pp.charAt(0));
            }else{
                if(pp.length() == 0){

                }else{
                    throw new RuntimeException("java implementation error.");
                }
            }
        }
    }
}
