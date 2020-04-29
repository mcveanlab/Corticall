package uk.ac.ox.well.cortexjdk;

public class CortexJDK {
    public static final String progName = "Corticall";
    public static final String progDesc = "germline DNM caller based on LdBG assembly of pathogens";
    public static final String rootPackage = "uk.ac.ox.well.cortexjdk";

    public static void main(String[] args) throws Exception {
        Main.start(progName, progDesc, rootPackage, args);
    }
}