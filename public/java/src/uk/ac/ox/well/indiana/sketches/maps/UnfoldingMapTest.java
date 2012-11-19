package uk.ac.ox.well.indiana.sketches.maps;

import com.modestmaps.StaticMap;
import com.modestmaps.core.Point2f;
import com.modestmaps.geo.Location;
import com.modestmaps.providers.Microsoft;
import processing.core.PImage;
import uk.ac.ox.well.indiana.sketches.Sketch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;

public class UnfoldingMapTest extends Sketch {
    @Argument(fullName="nonsense", shortName="n", doc="Nonsense argument")
    public String NONSENSE;

    /*
    private UnfoldingMap map;

    public void setup() {
        size(800, 600, GLConstants.GLGRAPHICS);
        map = new UnfoldingMap(this);
        MapUtils.createDefaultEventDispatcher(this, map);
    }

    public void draw() {
        map.draw();
    }
    */

    public void setup() {
        size(800, 600);
        noLoop();
    }

    public void draw() {
        StaticMap m = new StaticMap(this, new Microsoft.AerialProvider(), new Point2f(width/2, height), new Location(51.5f, -0.137f), 12);

        PImage img = m.draw(true);
        image(img,0,0);

        //img = atkinsonDither(img);
        image(img,width/2,0);

        println("done");

    }
}

