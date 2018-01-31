package uk.ac.ox.well.cortexjdk.utils.alignment.sw;

class SWCell{
    double score[] = new double[3];
    int prevType[] = new int[3];
    public void setScoreAt(int type,double sc){
        score[type] = sc;
    }
    public double getScoreAt(int type){
        return score[type];
    }
    public int getPrevTypeAt(int type){
        return prevType[type];
    }
    public void setPrevTypeAt(int ctype,int ptype){
        prevType[ctype] = ptype;
    }


}
