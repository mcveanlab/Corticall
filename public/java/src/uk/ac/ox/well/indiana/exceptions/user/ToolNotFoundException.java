package uk.ac.ox.well.indiana.exceptions.user;

public class ToolNotFoundException extends UserException {
    public ToolNotFoundException(String toolName) {
        super("The tool '" + toolName + "' was not found.");
    }

    public ToolNotFoundException(String toolName, Throwable e) {
        super("The tool '" + toolName + "' was not found: " + e.getMessage() + ".", e);
    }
}
