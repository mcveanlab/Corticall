package uk.ac.ox.well.cortexjdk;

public class CortexJDK {
    public static final String progName = "CortexJDK";
    public static final String progDesc = "tools for manipulating (Mc)Cortex de-novo assembly graph and link data";
    public static final String rootPackage = "uk.ac.ox.well.cortexjdk";

    public static void main(String[] args) throws Exception {
        Main.start(progName, progDesc, rootPackage, args);
    }
}