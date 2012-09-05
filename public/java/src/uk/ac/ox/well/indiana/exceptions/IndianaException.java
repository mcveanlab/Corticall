package uk.ac.ox.well.indiana.exceptions;

public class IndianaException extends Exception {
    public IndianaException(String msg) { super("\n        " + msg + "\n\nStack trace:"); }
    public IndianaException(String msg, Throwable e) { super("\n        " + msg + "\n\nStack trace:", e); }
    private IndianaException(Throwable e) { super("\n\nStack trace:", e); }
}
