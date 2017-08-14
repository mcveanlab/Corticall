package uk.ac.ox.well.cortexjdk.commands;

import org.slf4j.Logger;
import uk.ac.ox.well.cortexjdk.Main;
import uk.ac.ox.well.cortexjdk.utils.arguments.ArgumentHandler;

public abstract class Module implements Command {
    public Logger log = Main.getLogger();
    public String[] args;

    public Module() {}

    public void init() {
        ArgumentHandler.parse(this, args);
    }

    public abstract void execute();
}
