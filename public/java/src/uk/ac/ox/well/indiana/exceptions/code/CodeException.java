package uk.ac.ox.well.indiana.exceptions.code;

import uk.ac.ox.well.indiana.exceptions.IndianaException;

public class CodeException extends IndianaException {
    public CodeException(String msg) {
        super(msg);
    }

    public CodeException(String msg, Throwable e) {
        super(msg, e);
    }
}
