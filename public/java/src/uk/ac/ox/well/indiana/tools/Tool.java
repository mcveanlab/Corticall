package uk.ac.ox.well.indiana.tools;

import uk.ac.ox.well.indiana.IndianaModule;
import uk.ac.ox.well.indiana.utils.arguments.ArgumentParser;

public abstract class Tool implements IndianaModule {
    public String[] args;

    public Tool() {}

    public void init() {
        ArgumentParser.parse(this, args);
    }

    public abstract int execute();
}
