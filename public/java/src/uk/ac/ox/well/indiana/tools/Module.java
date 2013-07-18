package uk.ac.ox.well.indiana.tools;

import org.slf4j.Logger;
import uk.ac.ox.well.indiana.Indiana;
import uk.ac.ox.well.indiana.IndianaModule;
import uk.ac.ox.well.indiana.utils.arguments.ArgumentHandler;

public abstract class Module implements IndianaModule {
    public Logger log = Indiana.getLogger();
    public String[] args;

    public Module() {}

    public void init() {
        ArgumentHandler.parse(this, args);
    }

    public abstract void execute();
}
