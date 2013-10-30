package uk.ac.ox.well.indiana;

public class IndianaException extends RuntimeException {
    public IndianaException(final String message) {
        super(message);
    }

    public IndianaException(final String message, final Throwable throwable) {
        super(message, throwable);
    }
}
