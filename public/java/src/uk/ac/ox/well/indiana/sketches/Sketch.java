package uk.ac.ox.well.indiana.sketches;

import org.slf4j.Logger;
import processing.core.PApplet;
import uk.ac.ox.well.indiana.IndianaMain;
import uk.ac.ox.well.indiana.IndianaModule;
import uk.ac.ox.well.indiana.utils.arguments.ArgumentHandler;

public abstract class Sketch extends PApplet implements IndianaModule {
    public Logger log = IndianaMain.getLogger();

    public Sketch() {}

    public void init() {
        ArgumentHandler.parse(this, args);

        initialize();

        super.init();
    }

    public void initialize() {}
}
