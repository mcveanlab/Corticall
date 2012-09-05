package uk.ac.ox.well.indiana.exceptions.code;

public class ToolResolutionException extends CodeException {
    public ToolResolutionException(String toolName) {
        super("The tool '" + toolName + "' could not be instantiated.");
    }

    public ToolResolutionException(String toolName, Throwable e) {
        super("The tool '" + toolName + "' could not be instantiated: " + e.getMessage() + ".", e);
    }
}
