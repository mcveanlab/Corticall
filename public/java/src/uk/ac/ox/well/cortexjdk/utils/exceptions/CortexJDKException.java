package uk.ac.ox.well.cortexjdk.utils.exceptions;

public class CortexJDKException extends RuntimeException {
    public CortexJDKException(final String message) {
        super(message);
    }

    public CortexJDKException(final String message, final Throwable throwable) {
        super(message, throwable);
    }
}
