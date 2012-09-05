package uk.ac.ox.well.indiana.exceptions.user;

import uk.ac.ox.well.indiana.exceptions.IndianaException;

public class UserException extends IndianaException {
    public UserException(String msg) {
        super(msg);
    }

    public UserException(String msg, Throwable e) {
        super(msg, e);
    }
}
