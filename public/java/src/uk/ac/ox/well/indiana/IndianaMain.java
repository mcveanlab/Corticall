package uk.ac.ox.well.indiana;

import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.ClassLocator;

public class IndianaMain {
    public static void main(String[] args) {
        System.out.println("Hello World!");

        ClassLocator<Tool> locator = new ClassLocator<Tool>(Tool.class);
    }
}