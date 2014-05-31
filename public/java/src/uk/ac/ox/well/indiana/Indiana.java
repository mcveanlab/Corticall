package uk.ac.ox.well.indiana;

public class Indiana {
    public static final String progName = "INDIANA";
    public static final String progDesc = "tools for constructing and manipulating de-novo assembly data";
    public static final String rootPackage = "uk.ac.ox.well.indiana";
    public static final String commandPackage = "uk.ac.ox.well.indiana.commands";

    public static void main(String[] args) throws Exception {
        Main.start(progName, progDesc, rootPackage, commandPackage, args);
    }
}