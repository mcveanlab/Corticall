package uk.ac.ox.well.indiana.commands;

import org.slf4j.Logger;
import processing.core.PApplet;
import uk.ac.ox.well.indiana.Indiana;
import uk.ac.ox.well.indiana.Main;
import uk.ac.ox.well.indiana.utils.arguments.ArgumentHandler;

public abstract class Sketch extends PApplet implements Command {
    public Logger log = Main.getLogger();

    public Sketch() {}

    public void init() {
        ArgumentHandler.parse(this, args);

        initialize();

        super.init();
    }

    public void initialize() {}
}
